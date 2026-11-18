#!/bin/bash -l
#Usage # sbatch launch_baseline_signature_knn_A549.sh [run ...]      (from scripts/predict_targets)
#   Signature-similarity baseline (baseline_signature_knn.py; STATUS 6.0b item 1): for each run given (default: the
#   current recipe, MetChange / expr / TF _qsigtop) and its _shuff control, the k most similar training drugs vote
#   their targets, scored with the model's own ranking code on the same folds, queries and candidates. Two inputs:
#     sources  exactly the model's input (selected features, signs, on network nodes): does the network add anything?
#     profile  the full continuous signature (Spearman): is there more in the data than the model's input keeps?
#   Then compare_runs.py for each baseline pair, and the performance figure with the baselines next to the model.
#   Must run here: each fold is checked against the run (held-out drugs and test perturbations rebuilt from the
#   run's own inputs) and the script stops if the inputs differ, which they do on the laptop (older labels, no
#   gene network).
#   Variants: K="10 25 50" sbatch ... (k sweep); REPS=sources sbatch ... (one representation).
#   Output: <runs_dir>/<run>__knn_<rep>_k<k>{,_shuff}/ (run-shaped: cv_drug_ranks*.tsv, fold<i>/, compare_real_vs_shuffled/)
#           <runs_dir>/performance_figure/vs_baseline/performance_panels_<metric>.{png,pdf,tsv}
#SBATCH --job-name=baseline_knn
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 4
#SBATCH --mem=48GB
#SBATCH --time=0-04:00:00
#SBATCH --output=./baseline_knn_%j.out
#SBATCH --error=./baseline_knn_%j.err
#SBATCH -p batch
#SBATCH --qos=normal

project_dir_scratch="/scratch/users/gsantamaria/test_large_files/metabognn_dev"
cell="A549"
runs_dir="${project_dir_scratch}/results/predict_targets/${cell}"

runs=("$@")
[[ ${#runs[@]} -gt 0 ]] || runs=(metchange_rank_noGlobals_lr1e-3_min15_qsigtop expr_noGlobals_lr1e-3_min15_qsigtop
                                 tf_active200_noGlobals_lr1e-3_min15_qsigtop)
read -r -a ks <<< "${K:-25}"
read -r -a reps <<< "${REPS:-sources profile}"

conda activate metabolic_drug_repurposing
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
status=0
done_pairs=()
for r in "${runs[@]}"; do
    for side in "" "_shuff"; do
        d="${runs_dir}/${r}${side}"
        if [[ ! -d "$d" ]]; then echo "WARNING: ${d} not found, skipped"; status=1; continue; fi
        for rep in "${reps[@]}"; do
            for k in "${ks[@]}"; do
                echo "=== ${r}${side}: ${rep}, k ${k}"
                python baseline_signature_knn.py --run-dir "$d" --representation "$rep" --k "$k" || status=1
            done
        done
    done
    for rep in "${reps[@]}"; do
        for k in "${ks[@]}"; do
            b="${runs_dir}/${r}__knn_${rep}_k${k}"
            if [[ -s "${b}/cv_drug_ranks.tsv" && -s "${b}_shuff/cv_drug_ranks.tsv" ]]; then
                python compare_runs.py --run real="$b" --run shuffled="${b}_shuff" --control shuffled \
                    --out-dir "${b}/compare_real_vs_shuffled" > "${b}/compare.log" 2>&1 || status=1
                done_pairs+=("${r}|${rep}|${k}")
            fi
        done
    done
done

# figure: the model rows plus, per signature, each baseline that finished (first k only, to keep it readable)
label() { case "$1" in metchange*) echo mc;; expr*) echo expr;; tf*) echo tf;; esac; }
name() { case "$1" in mc) echo MetChange;; expr) echo expression;; tf) echo TF;; esac; }
extra=()
for e in "${done_pairs[@]}"; do
    IFS="|" read -r r rep k <<< "$e"
    [[ "$k" == "${ks[0]}" ]] || continue
    s=$(label "$r")
    what=$([[ "$rep" == sources ]] && echo "kNN, same input" || echo "kNN, full profile")
    extra+=(--extra "${s}=${r}__knn_${rep}_k${k}:$(name "$s") ${what}")
done
if [[ ${#extra[@]} -gt 0 ]]; then
    for m in rank_probc_heldout spec_prob_timedose_heldout rank_prob_heldout; do
        python plot_performance_panels.py --runs-dir "$runs_dir" --out-dir "${runs_dir}/performance_figure/vs_baseline" \
            --metric "$m" "${extra[@]}" || status=1
    done
fi
conda deactivate
exit $status
