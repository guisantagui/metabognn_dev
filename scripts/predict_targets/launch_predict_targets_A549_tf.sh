#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=predict_targets_A549_tf
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH -c 16
#SBATCH --mem=96GB
#SBATCH --time=0-12:00:00
#Define sdout path
#SBATCH --output=./predict_targets_A549_tf.out
#Define sderr path
#SBATCH --error=./predict_targets_A549_tf.err
# GPU partition and request: adjust to the cluster (iris: -p gpu with -G 1)
#SBATCH -p gpu
#SBATCH -G 1
#SBATCH --qos=normal

# Python env with torch (CUDA), numpy, pandas, scipy, scikit-learn, networkx
#conda activate torch_gnn
conda activate metabolic_drug_repurposing

project_dir="/home/users/gsantamaria/projects/metabognn_dev"
project_dir_scratch="/scratch/users/gsantamaria/test_large_files/metabognn_dev/"
cell="A549"
net_dir="${project_dir}/data/networks/tf"
cell_dir="${project_dir_scratch}/data/perturbation/lincs/cell_matrices/${cell}"   # as in the other launchers
hgnc_f="${project_dir}/data/utilities/hgnc_complete_set.txt"
#labels_dir="${project_dir}/results/labels/${cell}/expr"   # limma_DE.R on the expression
labels_dir="${project_dir_scratch}/results/labels/${cell}/expr"
out_dir="${project_dir}/results/predict_targets/${cell}/tf"
device="gpu" # cpu, cuda

# Control baseline (baseline_ctrl.tsv) and perturbation table (perturbations.tsv), derived from expr.csv and
# mdat.csv: (re)computed when missing or older than their inputs
if [ ! -f "${cell_dir}/baseline_ctrl.tsv" ] || [ ! -f "${cell_dir}/perturbations.tsv" ] \
   || [ "${cell_dir}/expr.csv" -nt "${cell_dir}/baseline_ctrl.tsv" ] \
   || [ "${cell_dir}/mdat.csv" -nt "${cell_dir}/perturbations.tsv" ]; then
        python prepare_inputs.py --cell-dir "$cell_dir" --hgnc "$hgnc_f"
fi

# Ablation (Menacho et al. 2026 setting): differentially expressed TFs propagated on the protein network
# without transcriptional edges (build_network.py --network tf).
# Signature rows (gene symbols or Ensembl ids) are mapped to the network's HGNC symbols with --hgnc.
# Relation weights default to build_network.SUGGESTED_ALPHA; ablations: --alpha grn_binding=0,
# --alpha second_messenger=0, --no-vae.
# Cross-validation by drug: --n-folds 5 runs all folds in this job; or submit as an array (#SBATCH --array=0-4)
# with --n-folds 5 --fold "$SLURM_ARRAY_TASK_ID" (the job that completes the last fold pools them).
# --linear-only: logistic regression only, minutes instead of hours.
python metabolic_target_prediction.py \
        --nodes "${net_dir}/nodes.tsv" \
        --edges "${net_dir}/edges.tsv" \
        --baseline "${cell_dir}/baseline_ctrl.tsv" \
        --perturbations "${cell_dir}/perturbations.tsv" \
        --signatures "${labels_dir}/log2FC.csv" \
        --qvalues "${labels_dir}/qvals.csv" \
        --mode tf \
        --hgnc "$hgnc_f" \
        --drug-targets "${project_dir}/data/drug_targets/combined/drug_targets.tsv" \
        --device "$device" \
        --predict-unlabeled \
        --n-folds \
        --no-vae \
        --out-dir "${out_dir}"

conda deactivate
