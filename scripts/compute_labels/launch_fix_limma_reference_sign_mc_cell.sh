#!/bin/bash --login
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name="corr_sign_mc_A549"
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32GB
#SBATCH -c 1
#SBATCH --time=00-03:00:00
#Define sdout path
#SBATCH --output="logs/corr_sign_mc_A549.out"
#Define sderr path
#SBATCH --error="logs/corr_sign_mc_A549.err"
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

proj_dir="${1:-/scratch/users/gsantamaria/test_large_files/metabognn_dev/}"
cell="${2:-A549}"
penalty="${3:-rank}"

in_dir="${proj_dir}/results/labels/${cell}/metchange/${penalty}/"
out_dir="${proj_dir}/results/labels/${cell}/metchange/${penalty}_signfixed/"
mdat="${proj_dir}/data/perturbation/lincs/cell_matrices/${cell}/mdat.csv"

mkdir -p "$out_dir"

# R's own level order, in the same environment as the limma jobs (conda r-4.3.1, login shell, same locale):
# which perturbations limma used as the reference level, i.e. got control - treated
conda activate r-4.3.1
Rscript limma_reference_order.R "$mdat" "${out_dir}/flipped_R.tsv" || { echo "limma_reference_order.R failed"; exit 1; }
conda deactivate

export MAMBA_ROOT_PREFIX="/work/projects/car_t_methylome/software/micromamba"

eval "$(micromamba shell hook --shell bash)"

micromamba activate metabolic_drug_repurposing

# conda / micromamba activation scripts do not run cleanly under set -e, so each step checks its own exit status
python fix_limma_reference_sign.py --dir "$in_dir" --out-dir "$out_dir" --name-flips "${out_dir}/flipped_R.tsv" \
    || { echo "fix_limma_reference_sign.py failed"; exit 1; }

micromamba deactivate