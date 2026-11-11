#!/bin/bash --login
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=mc_concat_celllist
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64GB
#SBATCH -c 1
#SBATCH --time=00-01:00:00
#Define sdout path
#SBATCH --output=./logs/mc_concat_celllist.out
#Define sderr path
#SBATCH --error=./logs/mc_concat_celllist.err
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

set -euo pipefail

micromamba activate metabolic_drug_repurposing

cell_file="/home/users/gsantamaria/projects/metabognn_dev/data/LINCS/core_cell_lines.tsv"
project_dir="/home/users/gsantamaria/projects/metabognn_dev/"
project_dir_scratch="/scratch/users/gsantamaria/test_large_files/metabognn_dev/"
penalty="rank"

mc_dir="${project_dir_scratch}/results/metChange/"
overwrite=false

while IFS=$'\t' read -r cell; do
    mc_dir_cell="${mc_dir}/${cell}/${penalty}/"
    out_dir_cell="${mc_dir}/${cell}/${penalty}_concat/"
    mkdir -p "$out_dir_cell"
    echo "Concatenating ${cell} MetChange results..."
    python concat_mc_res.py --mc_dir "$mc_dir_cell" --out_dir "$out_dir_cell"
done < "$cell_file"

micromamba deactivate