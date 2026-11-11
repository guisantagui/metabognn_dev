#!/bin/bash -l
#Usage # sbatch launch_compare_runs.sh --run real=<dir> --run shuffled=<dir> [--control shuffled] [...]
#   (all arguments go to compare_runs.py; submitted by submit_with_shuffle_control.sh after both runs end)
#SBATCH --job-name=compare_runs
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 2
#SBATCH --mem=8GB
#SBATCH --time=0-00:30:00
#SBATCH --output=./compare_runs_%j.out
#SBATCH --error=./compare_runs_%j.err
#SBATCH -p batch
#SBATCH --qos=normal

conda activate metabolic_drug_repurposing
python compare_runs.py "$@"
conda deactivate
