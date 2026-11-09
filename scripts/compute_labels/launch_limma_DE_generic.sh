#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=200GB
#SBATCH -c 128
#SBATCH --time=02-00:00:00
#Define sdout path
#Define sderr path
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal

conda activate r-4.3.1
# One BLAS thread per worker: the workers already use all the cores
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# Results are checkpointed in <out_dir>/checkpoints: if the job hits the time
# limit, resubmit it as is and it continues where it stopped.

# Script arguments
################################################################################

project_dir="$1"
cell="$2"
workers="$3"


# Get paths
################################################################################

echo "Running limma of ${cell}."

cell_mat_dir="${project_dir}/data/perturbation/lincs/cell_matrices/"
out_dir="${project_dir}/results/labels/"

mkdir -p "$out_dir"

expr="${cell_mat_dir}/${cell}/expr_ensembl.csv"
mdat="${cell_mat_dir}/${cell}/mdat.csv"
out_dir_cell="${out_dir}/${cell}/expr/"

mkdir -p "$out_dir_cell"

# Run limma
################################################################################

Rscript limma_DE.R --expr "$expr" \
        --mdat "$mdat" \
        --workers "$workers" \
        --out_dir "$out_dir_cell"
conda deactivate