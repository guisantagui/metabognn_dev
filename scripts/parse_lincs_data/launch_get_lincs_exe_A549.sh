#!/bin/bash
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=lincs_A549
#SBATCH --mem=150G
#SBATCH --time=00-15:00:00
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#Define sdout path
#SBATCH --output=./launch_get_lincs_A549.out
#Define sderr path
#SBATCH --error=./launch_get_lincs_A549.err
#Define the queue (Quality Of Service) to which the task shall be submitted to
#SBATCH --partition=batch

gctxFile_pert="/scratch/users/gsantamaria/test_large_files/data/perturbation/lincs/cp_predicted_RNAseq_profiles.gctx"
gctxFile_ctrl="/scratch/users/gsantamaria/test_large_files/data/perturbation/lincs/ctl_predicted_RNAseq_profiles.gctx"
cell="A549"
ctrls_keep="DMSO"
hgnc="/home/users/gsantamaria/projects/metabo_pert/data/utilities/hgnc_complete_set.txt"
out_dir="/scratch/users/gsantamaria/test_large_files/data/perturbation/lincs/cell_matrices/"

out_dir_cell="${out_dir}/${cell}/"

mkdir -p "$out_dir_cell"


conda activate r-4.3.1

Rscript get_lincs_exe.R \
    --gctx_pert "$gctxFile_pert" \
    --gctx_ctrl "$gctxFile_ctrl" \
    --cell "$cell" \
    --ctrls_keep "$ctrls_keep" \
    --out_dir "$out_dir_cell"

conda deactivate

# Convert symbols to ensembl_gene_id
export MAMBA_ROOT_PREFIX="/work/projects/car_t_methylome/software/micromamba"
eval "$(micromamba shell hook --shell bash)"
micromamba activate metabolic_drug_repurposing

expr_f = args.expr
hgnc_f = args.hgnc
symbol_col = args.symbol_col
is_log2 = args.is_log2
delete_original = args.delete_original
out_dir = args.out_dir

python expr_symb_to_ensembl.py \
    --expr "${out_dir_cell}expr.csv" \
    --hgnc "$hgnc" \
    --symbol_col "V1" \
    --is_log2 \
    --delete_original \
    --out_dir "$out_dir_cell"
micromamba deactivate