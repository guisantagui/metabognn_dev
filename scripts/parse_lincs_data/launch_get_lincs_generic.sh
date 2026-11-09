#!/bin/bash
#Usage # sbatch [this script]
#Name of the job
#SBATCH --mem=150G
#SBATCH --time=00-15:00:00
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#Define sdout path
#Define sderr path
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH --partition=batch

# Input arguments
################################################################################

pert_file="$1"
ctrl_file="$2"
cell="$3"
ctrls_keep="$4"
hgnc="$5"
gem_f="$6"
proj_dir="$7"


# Setting up output directory
################################################################################

out_dir="${proj_dir}/data/perturbation/lincs/cell_matrices/"

out_dir_cell="${out_dir}/${cell}/"

mkdir -p "$out_dir_cell"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate r-4.3.1

# Extract data from LINCS
################################################################################

Rscript get_lincs_exe.R \
    --gctx_pert "$pert_file" \
    --gctx_ctrl "$ctrl_file" \
    --cell "$cell" \
    --ctrls_keep "$ctrls_keep" \
    --out_dir "$out_dir_cell"

conda deactivate

# Translate symbols to ENSEMBL
################################################################################

# Convert symbols to ensembl_gene_id
export MAMBA_ROOT_PREFIX="/work/projects/car_t_methylome/software/micromamba"
eval "$(micromamba shell hook --shell bash)"
micromamba activate metabolic_drug_repurposing

python expr_symb_to_ensembl.py \
    --expr "${out_dir_cell}expr.csv" \
    --hgnc "$hgnc" \
    --gem "$gem_f" \
    --symbol_col "V1" \
    --is_log2 \
    --out_dir "$out_dir_cell"
micromamba deactivate