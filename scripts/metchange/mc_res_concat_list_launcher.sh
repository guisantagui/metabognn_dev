#!/bin/bash --login
#Usage # sbatch [this script]
# One array task per cell line of cell_file (line N = task N): set --array to 1-<number of lines>.
#Name of the job
#SBATCH --job-name=mc_concat_celllist
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=64GB
#SBATCH -c 8
#SBATCH --time=00-01:00:00
#SBATCH --array=1-9
#Define sdout path
#SBATCH --output=./logs/mc_concat_celllist_%a.out
#Define sderr path
#SBATCH --error=./logs/mc_concat_celllist_%a.err
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

# polars reads each CSV with all the cores of the task
export POLARS_MAX_THREADS="${SLURM_CPUS_PER_TASK:-1}"

cell=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$cell_file" | cut -f1)
if [ -z "$cell" ]; then
    echo "No cell line on line ${SLURM_ARRAY_TASK_ID} of ${cell_file}"
    exit 1
fi

mc_dir_cell="${mc_dir}/${cell}/${penalty}/"
out_dir_cell="${mc_dir}/${cell}/${penalty}_concat/"
if [ "$overwrite" = false ] && ls "${out_dir_cell}"/*_concat.csv >/dev/null 2>&1; then
    echo "${cell}: already concatenated in ${out_dir_cell} (set overwrite=true to redo)"
    exit 0
fi
mkdir -p "$out_dir_cell"
echo "Concatenating ${cell} MetChange results..."
python concat_mc_res.py --mc_dir "$mc_dir_cell" --out_dir "$out_dir_cell"

micromamba deactivate
