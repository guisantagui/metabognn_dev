#!/bin/bash --login
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=delta_rnkAll
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32GB
#SBATCH -c 1
#SBATCH --time=00-01:00:00
#Define sdout path
#SBATCH --output=delta_A549_rnkAll.out
#Define sderr path
#SBATCH --error=delta_A549_rnkAll.err
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

# Concatenates all the metchange files and computes deltas and averages, generating corresponging
# metadata files

eval "$(micromamba shell hook --shell bash)"
micromamba activate metabolic_drug_repurposing

in_dir="/scratch/users/gsantamaria/test_large_files/metabo_perts/results/metChange/"
cell_mat_dir=/scratch/users/gsantamaria/test_large_files/data/perturbation/lincs/cell_matrices/
cell="A549"
penalty="rank"
batch_col="batch"
out_dir="/scratch/users/gsantamaria/test_large_files/metabo_perts/results/mc_deltas/"


mc_dir="${in_dir}/${cell}/${penalty}/"
mdat="${cell_mat_dir}/${cell}/mdat.csv"
out_dir_cell="${out_dir}/${cell}/${penalty}/"

mkdir -p "$out_dir_cell"

python get_delta_mat.py \
    --mc_dir "$mc_dir" \
    --mdat "$mdat" \
    --batch_col "$batch_col" \
    --out_dir "$out_dir_cell"
micromamba deactivate