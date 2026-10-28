#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=limma_A549_expr
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=128GB
#SBATCH -c 128
#SBATCH --time=01-00:00:00
#Define sdout path
#SBATCH --output=./limma_A549_expr.out
#Define sderr path
#SBATCH --error=./limma_A549_expr.err
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

conda activate r-4.3.1

cell="A549"
cell_mat_dir=/scratch/users/gsantamaria/test_large_files/data/perturbation/lincs/cell_matrices/
workers=64
out_dir="/home/users/gsantamaria/projects/metabo_pert/results/labels/"

mkdir -p "$out_dir"

expr="${cell_mat_dir}/${cell}/expr_ensembl.csv"
mdat="${cell_mat_dir}/${cell}/mdat.csv"
out_dir_cell="${out_dir}/${cell}/expr/${penalty}/"

mkdir -p "$out_dir_cell"

Rscript limma_DE.R --expr "$expr" \
        --mdat "$mdat" \
        --workers "$workers" \
        --out_dir "$out_dir_cell"
conda deactivate