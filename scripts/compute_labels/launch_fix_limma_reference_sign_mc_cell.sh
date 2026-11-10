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

export MAMBA_ROOT_PREFIX="/work/projects/car_t_methylome/software/micromamba"

eval "$(micromamba shell hook --shell bash)"

micromamba activate metabolic_drug_repurposing

#proj_dir="$1"
#cell="$2"

proj_dir="/scratch/users/gsantamaria/test_large_files/metabognn_dev/"
cell="A549"
penalty="rank"

in_dir="${proj_dir}/results/labels/${cell}/metchange/${penalty}/"
out_dir="${proj_dir}/results/labels/${cell}/metchange/${penalty}_signfixed/"

mkdir -p "$out_dir"

python fix_limma_reference_sign.py --dir "$in_dir" --out-dir "$out_dir"

micromamba deactivate