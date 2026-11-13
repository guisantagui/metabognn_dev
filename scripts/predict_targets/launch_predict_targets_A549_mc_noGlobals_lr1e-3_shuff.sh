#!/bin/bash -l
#Usage # sbatch [this script]
#Name of the job
#SBATCH --job-name=predict_targets_A549_mc_noGlobals_lr1e-3_shuff
#SBATCH -N 1
#SBATCH --mail-user=guillem.santamaria@uni.lu
#SBATCH --mail-type=begin,end,fail
#SBATCH --ntasks-per-node=1
#SBATCH -c 16
#SBATCH --mem=96GB
#SBATCH --time=0-04:00:00
# One fold per array task (~30 min per fold at lr 5e-3; up to 2x the epochs here): short requests start sooner
#SBATCH --array=0-4
#Define sdout path
#SBATCH --output=./predict_targets_A549_mc_noGlobals_lr1e-3_shuff_%a.out
#Define sderr path
#SBATCH --error=./predict_targets_A549_mc_noGlobals_lr1e-3_shuff_%a.err
# GPU partition and request: adjust to the cluster (iris: -p gpu with -G 1)
#SBATCH -p gpu
#SBATCH -G 1
#SBATCH --qos=normal

# Python env with torch (CUDA), numpy, pandas, scipy, scikit-learn, networkx
#conda activate torch_gnn
conda activate metabolic_drug_repurposing

project_dir="/home/users/gsantamaria/projects/metabognn_dev"
project_dir_scratch="/scratch/users/gsantamaria/test_large_files/metabognn_dev"
cell="A549"
penalty="rank"
net_dir="${project_dir}/data/networks/combined"
cell_dir="${project_dir_scratch}/data/perturbation/lincs/cell_matrices/${cell}"   # as in the other launchers
hgnc_f="${project_dir}/data/utilities/hgnc_complete_set.txt"
labels_dir="${project_dir}/results/labels/${cell}/metchange/${penalty}"   # rerun limma_DE.R (fixed) first
labels_dir_expr="${project_dir_scratch}/results/labels/${cell}/expr"
out_dir="${project_dir}/results/predict_targets/${cell}/metchange_${penalty}_noGlobals_lr1e-3_shuff"
device="cuda" # cpu, cuda

# Control baseline (baseline_ctrl.tsv) and perturbation table (perturbations.tsv), derived from expr.csv and
# mdat.csv: (re)computed when missing or older than their inputs
if [ ! -f "${cell_dir}/baseline_ctrl.tsv" ] || [ ! -f "${cell_dir}/perturbations.tsv" ] \
   || [ "${cell_dir}/expr.csv" -nt "${cell_dir}/baseline_ctrl.tsv" ] \
   || [ "${cell_dir}/mdat.csv" -nt "${cell_dir}/perturbations.tsv" ]; then
        python prepare_inputs.py --cell-dir "$cell_dir" --hgnc "$hgnc_f"
fi

# MetChange inconsistency is sign-flipped so that positive = more production (Schuhknecht benchmark convention).
# Relation weights default to build_network.SUGGESTED_ALPHA; ablations: --alpha grn_binding=0,
# --alpha second_messenger=0, --no-vae.
# Cross-validation by drug: --n-folds 5 runs all folds in this job; or submit as an array (#SBATCH --array=0-4)
# with --n-folds 5 --fold "$SLURM_ARRAY_TASK_ID" (the job that completes the last fold pools them).
# --linear-only: logistic regression only, minutes instead of hours.
# Array job: task i runs fold i (same split and seed as running all folds in one job); the task that finishes
# last pools the folds (cv_*). If cv_* is missing (e.g. a task failed and was resubmitted alone with
# sbatch --array=<i>), pool with --aggregate-only --n-folds 5 --out-dir <out_dir>. compare_runs.py reads the fold
# files and does not need cv_*.
# --lr 0.001 (default 0.005), --patience 10 (default 5), --ffn-epochs 60 (default 30): training stopped after
# 1-3 epochs in several folds of metchange_rank_noGlobals, and some members kept epoch 26 of 30 in other runs
# (reports/2026-11-13_A549_matched_negatives.md, 5.5).
python metabolic_target_prediction.py \
        --nodes "${net_dir}/nodes.tsv" \
        --edges "${net_dir}/edges.tsv" \
        --baseline "${cell_dir}/baseline_ctrl.tsv" \
        --perturbations "${cell_dir}/perturbations.tsv" \
        --signatures "${labels_dir}/log2FC.csv" \
        --qvalues "${labels_dir}/qvals.csv" \
        --signature-sign -1 \
        --mode metchange \
        --drug-targets "${project_dir}/data/drug_targets/combined/drug_targets.tsv" \
        --drug-metabolites "${net_dir}/drug_metabolite_matches.tsv" \
        --device "$device" \
        --predict-unlabeled \
        --de-counts "${labels_dir_expr}/fit_info_nsig.csv" \
        --min-de-genes 200 \
        --score-active-only \
        --ignore-qvalue --top-n 100 \
        --n-folds 5 --fold "$SLURM_ARRAY_TASK_ID" \
        --lr 0.001 --patience 10 --ffn-epochs 60 \
        --shuffle-signatures \
        --no-vae \
        --no-globals \
        --out-dir "${out_dir}"

conda deactivate
