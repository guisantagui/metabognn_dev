#!/bin/bash --login
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=get_lincs_master
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16GB
#SBATCH -c 1
#SBATCH --time=00-01:00:00
#Define sdout path
#SBATCH --output=./logs/get_lincs_master.out
#Define sderr path
#SBATCH --error=./logs/get_lincs_master.err
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH -p batch
#SBATCH --qos=normal


cell_file="/home/users/gsantamaria/projects/metabognn_dev/data/LINCS/core_cell_lines.tsv"
pert_file="/scratch/users/gsantamaria/test_large_files/metabognn_dev/data/perturbation/lincs/cp_predicted_RNAseq_profiles.gctx"
ctrl_file="/scratch/users/gsantamaria/test_large_files/metabognn_dev/data/perturbation/lincs/ctl_predicted_RNAseq_profiles.gctx"
ctrls_keep="DMSO"
hgnc="/home/users/gsantamaria/projects/metabo_pert/data/utilities/hgnc_complete_set.txt"
gem_f="/home/users/gsantamaria/projects/metabo_pert/data/GEMs/Human2/Human-GEM/model/Human-GEM.xml"
proj_dir="/scratch/users/gsantamaria/test_large_files/metabognn_dev/"
overwrite=false
run_limma=true
limma_workers=64

mkdir -p "logs/"

while IFS=$'\t' read -r cell; do
    log2FC="${proj_dir}/results/labels/${cell}/expr/log2FC.csv"
    qvals="${proj_dir}/results/labels/${cell}/expr/qvals.csv"
    basal="${proj_dir}/results/labels/${cell}/expr/basal.csv"
    if [[ -f "$basal" && -f "$qvals" && -f "$log2FC" && $overwrite == false ]]; then
        echo "DE files already exist for ${cell}. Skipping."
    else
        echo "Launching ${cell} extraction..."
        sbatch \
            --job-name="lincs_${cell}" \
            --output="logs/get_lincs_${cell}.out" \
            --error="logs/get_lincs_${cell}.err" \
            launch_get_lincs_generic.sh \
            "$pert_file" \
            "$ctrl_file" \
            "$cell" \
            "$ctrls_keep" \
            "$hgnc" \
            "$gem_f" \
            "$proj_dir" \
            "$overwrite" \
            "$run_limma" \
            "$limma_workers"
    fi
done < "$cell_file"