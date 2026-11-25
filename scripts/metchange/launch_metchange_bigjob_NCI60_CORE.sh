#!/bin/bash --login
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=MC_master_NCI60_CORE
#SBATCH -N 1
#SBATCH --cpus-per-task=128
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=100GB
#SBATCH -c 1
#SBATCH --time=00-08:00:00
#Define sdout path
#SBATCH --output=MC_master_NCI60_CORE.out
#Define sderr path
#SBATCH --error=MC_master_NCI60_CORE.err
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

project_dir="/home/users/gsantamaria/projects/metabognn_dev/"
project_dir_scratch="/scratch/users/gsantamaria/test_large_files/metabognn_dev/"
cell="NCI60"
penalty="rank"
expr_f="/home/users/gsantamaria/projects/metabognn_dev/data/validation/human1_nci60/nci60_expr_ensembl.csv"
#cell_mat_dir="${project_dir_scratch}/data/perturbation/lincs/cell_matrices/"
#hgnc_f="${project_dir}/data/utilities/hgnc_complete_set.txt"
input_type="penalty"
batch_size=1000
gem_f="${project_dir}/data/GEMs/Human-GEM/model/Human-GEM.xml"
medium_f="${project_dir}/data/GEMs/gem_mediums/Human2-GEM/rpmi_medium.json"
out_dir="${project_dir_scratch}/results/mc_phenotype_assoc/NCI60_CORE/metchange/"

mkdir -p "$out_dir"


#expr_f="${cell_mat_dir}/${cell}/expr_ensembl.csv"
#mdat_f="${cell_mat_dir}/${cell}/mdat.csv"
out_dir_expr="$(dirname $out_dir)/${penalty}/"
out_dir_metchange="${out_dir}/${penalty}/"

mkdir -p "$out_dir_expr"
mkdir -p "$out_dir_metchange"

expr_pen_f="${expr_f%.csv}_rank.csv"
pen_base=$(basename "$expr_pen_f" .csv)
#mdat_base=$(basename "$mdat_f" .csv)

# Rank-based penalties: (re)computed only when missing or older than the expression
if [ ! -f "$expr_pen_f" ] || [ "$expr_f" -nt "$expr_pen_f" ]; then
    python create_penalties.py \
        --expr "$expr_f" \
        --threshold 2 \
        --pen_mode "$penalty" \
        --gem "$gem_f" \
        --out_dir $(dirname "$expr_f")
else
    echo "${expr_pen_f} is up to date, skipping create_penalties.py."
fi


python run_metchange_exe.py \
    --expr "$expr_pen_f" \
    --input_type "$input_type" \
    --gem "$gem_f" \
    --medium "$medium_f" \
    --n_workers 124 \
    --out_dir "$out_dir_metchange"