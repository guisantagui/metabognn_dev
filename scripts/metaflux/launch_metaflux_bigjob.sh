#!/bin/bash --login
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=MF_master_A549
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=200GB
#SBATCH -c 1
#SBATCH --time=00-08:00:00
#Define sdout path
#SBATCH --output=mf_master_A549.out
#Define sderr path
#SBATCH --error=mf_master_A549.err
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

# METAFlux on every LINCS sample of one cell line: GEM-gene batches of the
# expression (prepare_input_matrices.py, expression mode), then one array task
# per batch (launch_metaflux_array.sh). Afterwards, as for MetChange:
#   get_delta_mat.py --mc_dir "$out_dir_metaflux" ...   (concatenation, deltas)
#   limma_DE.R --expr <..._metaflux_<medium>_concat.csv> ...
# The environment needs clarabel (pip install clarabel), or osqp / gurobipy
# for --solver osqp / gurobi.

set -euo pipefail

module --force purge || true

module load env/development/2024a

export MAMBA_ROOT_PREFIX="/work/projects/car_t_methylome/software/micromamba"

eval "$(micromamba shell hook --shell bash)"

micromamba activate metabolic_drug_repurposing

cell="A549"
cell_mat_dir=/scratch/users/gsantamaria/test_large_files/data/perturbation/lincs/cell_matrices/
hgnc_f="/home/users/gsantamaria/projects/metabo_pert/data/utilities/hgnc_complete_set.txt"
batch_size=2000
gem_f="/home/users/gsantamaria/projects/metabo_pert/data/GEMs/Human2/Human-GEM/model/Human-GEM.xml"
medium_f="/home/users/gsantamaria/projects/metabo_pert/data/GEMs/gem_mediums/Human2-GEM/rpmi_medium.json"
# Quantile of each sample's positive MRAS given to reactions whose genes read 0 (0 = as the R package).
# With LINCS A549 some biomass-essential reactions read 0 in every sample: with 0 every flux is 0, and
# below ~0.25 the floor itself sets the biomass flux (see run_metaflux.py).
mras_floor=0.25
out_dir="/scratch/users/gsantamaria/test_large_files/metabo_perts/results/metaflux/"

expr_f="${cell_mat_dir}/${cell}/expr.csv"          # gene symbols (V1) x samples, log2
mdat_f="${cell_mat_dir}/${cell}/mdat.csv"
out_dir_expr="${cell_mat_dir}/batch_split/expression/"
out_dir_metaflux="${out_dir}/${cell}/floor_${mras_floor}/"
status_dir="${out_dir}/${cell}/floor_${mras_floor}_status/"

mkdir -p "$out_dir_expr" "$out_dir_metaflux" "$status_dir"

expr_base=$(basename "$expr_f" .csv)

# Remove batches from earlier runs, so that the array size below counts only these
rm -f "${out_dir_expr}/${expr_base}"_batch_*.csv "${out_dir_expr}"/$(basename "$mdat_f" .csv)_batch_*.csv

# Symbols -> ENSEMBL, GEM genes only, batches of samples (written in linear scale)
python ../metchange/prepare_input_matrices.py --expr "$expr_f" \
    --mdat "$mdat_f" \
    --input_type expression \
    --hgnc "$hgnc_f" \
    --batch_size "$batch_size" \
    --gem "$gem_f" \
    --out_dir "$out_dir_expr"

n_batches=$(ls "${out_dir_expr}/${expr_base}"_batch_*.csv | wc -l)
echo "Submitting ${n_batches} METAFlux batches"
mkdir -p logs

sbatch \
    --job-name="${SLURM_JOB_NAME}_arr" \
    --array=0-$((n_batches - 1))%8 \
    launch_metaflux_array.sh \
    "$gem_f" \
    "$medium_f" \
    "$out_dir_expr" \
    "$out_dir_metaflux" \
    "$expr_f" \
    "$mras_floor" \
    "$status_dir"
