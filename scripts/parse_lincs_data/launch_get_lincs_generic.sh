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
overwrite="$8"
run_limma="$9"
limma_workers="${10}"


# Setting up output directory
################################################################################

out_dir="${proj_dir}/data/perturbation/lincs/cell_matrices/"

out_dir_cell="${out_dir}/${cell}/"

mkdir -p "$out_dir_cell"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate r-4.3.1

# Extract data from LINCS
################################################################################

if [[ ! -f "${out_dir_cell}/expr.csv" || $overwrite == true ]]; then
    echo "Extracting ${cell} expression data from LINCS GCTX files..."
    Rscript get_lincs_exe.R \
        --gctx_pert "$pert_file" \
        --gctx_ctrl "$ctrl_file" \
        --cell "$cell" \
        --ctrls_keep "$ctrls_keep" \
        --out_dir "$out_dir_cell" \
        || { echo "LINCS extraction failed for ${cell}"; exit 1; }
else
    echo "expr.csv already exists in ${out_dir_cell}. Skipping extraction of ${cell} data."
fi

conda deactivate

# Translate symbols to ENSEMBL
################################################################################

# Convert symbols to ensembl_gene_id
export MAMBA_ROOT_PREFIX="/work/projects/car_t_methylome/software/micromamba"
eval "$(micromamba shell hook --shell bash)"

if [[ ! -f "${out_dir_cell}/expr_ensembl.csv" || $overwrite == true ]]; then
    echo "Converting ${cell} expr.csv HGNC gene symbols to ENSEMBL gene IDs..."
    micromamba activate metabolic_drug_repurposing

    python expr_symb_to_ensembl.py \
        --expr "${out_dir_cell}expr.csv" \
        --hgnc "$hgnc" \
        --gem "$gem_f" \
        --symbol_col "V1" \
        --is_log2 \
        --out_dir "$out_dir_cell" \
        || { echo "ENSEMBL conversion failed for ${cell}"; exit 1; }
    micromamba deactivate
else
    echo "expr_ensembl.csv already exists in ${out_dir_cell}. Skipping conversion to ENSEMBL gene IDs of ${cell} data..."
fi

# Compute DEGs with limma
################################################################################
log2FC="${proj_dir}/results/labels/${cell}/expr/log2FC.csv"
qvals="${proj_dir}/results/labels/${cell}/expr/qvals.csv"
basal="${proj_dir}/results/labels/${cell}/expr/basal.csv"
if [[ "$run_limma" == true ]] && \
   [[ "$overwrite" == true || ! -f "$basal" || ! -f "$qvals" || ! -f "$log2FC" ]]; then
    cd ../compute_labels

    mkdir -p "logs/"

    sbatch \
        --job-name="limma_${cell}_expr" \
        --output="logs/limma_${cell}_expr.out" \
        --error="logs/limma_${cell}_expr.err" \
        launch_limma_DE_generic.sh \
        "$proj_dir" \
        "$cell" \
        "$limma_workers"
else
    echo "Limma outputs already exist for ${cell}, or limma execution is disabled. Skipping."
fi