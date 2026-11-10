#!/bin/bash --login
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=targ_sim_A549
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=32GB
#SBATCH -c 1
#SBATCH --time=00-08:00:00
#Define sdout path
#SBATCH --output=mtarg_sim_A549.out
#Define sderr path
#SBATCH --error=targ_sim_A549.err
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal


export MAMBA_ROOT_PREFIX="/work/projects/car_t_methylome/software/micromamba"

eval "$(micromamba shell hook --shell bash)"

micromamba activate metabolic_drug_repurposing


project_dir="/home/users/gsantamaria/projects/metabognn_dev"
project_dir_scratch="/scratch/users/gsantamaria/test_large_files/metabognn_dev"
cell="A549"
cell_dir="${project_dir_scratch}/data/perturbation/lincs/cell_matrices/${cell}"   # as in the other launchers
labels_dir="${project_dir_scratch}/results/labels/${cell}/expr"
out_dir="${project_dir_scratch}/results/"
out_dir_cell="${out_dir}/signature_similarity/${cell}/"
min_de_genes=200

mkdir -p $out_dir_cell

python check_signature_target_similarity.py \
    --signatures ${labels_dir}/log2FC.csv \
    --perturbations ${cell_dir}/perturbations.tsv \
    --drug-targets ${project_dir}/data/drug_targets/combined/drug_targets.tsv \
    --de-counts ${labels_dir}/fit_info_nsig.csv --min-de-genes "$min_de_genes" \
    --out-dir "$out_dir_cell"