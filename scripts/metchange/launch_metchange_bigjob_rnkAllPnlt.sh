#!/bin/bash --login
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=MC_master_A549_rankAll
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=200GB
#SBATCH -c 1
#SBATCH --time=00-02:00:00
#Define sdout path
#SBATCH --output=mc_master_A549_rankAll.out
#Define sderr path
#SBATCH --error=mc_master_A549_rankAll.err
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

set -euo pipefail

module --force purge || true

module load env/development/2024a
module load math/Gurobi/12.0.1-GCCcore-13.3.0

export MAMBA_ROOT_PREFIX="/work/projects/car_t_methylome/software/micromamba"

eval "$(micromamba shell hook --shell bash)"

micromamba activate metabolic_drug_repurposing

cell="A549"
penalty="rank"
cell_mat_dir=/scratch/users/gsantamaria/test_large_files/data/perturbation/lincs/cell_matrices/
hgnc_f="/home/users/gsantamaria/projects/metabo_pert/data/utilities/hgnc_complete_set.txt"
input_type="penalty"
batch_size=1000
gem_f="/home/users/gsantamaria/projects/metabo_pert/data/GEMs/Human2/Human-GEM/model/Human-GEM.xml"
medium_f="/home/users/gsantamaria/projects/metabo_pert/data/GEMs/gem_mediums/Human2-GEM/rpmi_medium.json"
out_dir="/scratch/users/gsantamaria/test_large_files/metabo_perts/results/metChange/"

expr_f="${cell_mat_dir}/${cell}/expr_ensembl.csv"
mdat_f="${cell_mat_dir}/${cell}/mdat.csv"
out_dir_expr="${cell_mat_dir}/batch_split/${penalty}/"
out_dir_metchange="${out_dir}/${cell}/${penalty}/"

mkdir -p "$out_dir_expr"
mkdir -p "$out_dir_metchange"

# Create rank-based penalties
python create_penalties.py \
    --expr "$expr_f" \
    --hgnc "$hgnc_f" \
    --threshold 2 \
    --pen_mode "$penalty" \
    --gem "$gem" \
    --out_dir $(dirname "$expr_f")

expr_pen_f="${expr_f%.csv}_rank.csv"


python prepare_input_matrices.py --expr "$expr_pen_f" \
    --mdat "$mdat_f" \
    --input_type "$input_type" \
    --hgnc "$hgnc_f" \
    --batch_size "$batch_size" \
    --gem "$gem_f" \
    --out_dir "$out_dir_expr"

if [ $? -ne 0 ]; then
    echo "ERROR: prepare_input_matrices.py failed."
    exit 1
fi

max_fluxes_file="${out_dir_metchange}/lp1_fluxes.npy"

if [ ! -f "$max_fluxes_file" ]; then
    python run_metchange_exe_lp1.py --medium "$medium_f" \
        --gem "$gem_f" \
        --out_dir "$out_dir_metchange"
else
    echo "${max_fluxes_file} already exists."
fi

# ---------------------------------------------------------
# Submit MetChange array
# ---------------------------------------------------------

sbatch \
    --job-name="${SLURM_JOB_NAME}_arr" \
    --array=0-31%8 \
    launch_metchange_array.sh \
    "$gem_f" \
    "$medium_f" \
    "$out_dir_expr" \
    "$out_dir_metchange" \
    "$expr_f" \
    "$input_type"