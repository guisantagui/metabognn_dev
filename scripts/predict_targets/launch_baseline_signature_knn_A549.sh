#!/bin/bash -l
#Usage # sbatch launch_baseline_signature_knn_A549.sh                     (from scripts/predict_targets)
#       sbatch --array=0-<n-1> launch_baseline_signature_knn_A549.sh run1 ... run<n>   (other runs)
#   Signature-similarity baseline (baseline_signature_knn.py; STATUS 6.0b item 1). Array job, one task per run
#   (default: the current recipe, MetChange / expr / TF _qsigtop): each task scores the run AND its _shuff control,
#   both representations and every k in ONE call, reading the label files once:
#     sources  exactly the model's input (selected features, signs, on network nodes): does the network add anything?
#     profile  the full continuous signature (Spearman): is there more in the data than the model's input keeps?
#   then compare_runs.py for each baseline pair, then the performance figure with the baselines next to the model
#   (redrawn by every task with whatever pairs exist, so the last task to finish leaves the complete figure).
#   Must run here: each fold is checked against the run (held-out drugs and test perturbations rebuilt from the
#   run's own inputs) and the script stops if the inputs differ, which they do on the laptop.
#   Wall time is one label read per task: expr is the longest (~30-50 min), MetChange ~10-15, TF a few minutes.
#   Variants: K="10 25 50" sbatch ... (k sweep); REPS=sources sbatch ... (one representation).
#   Output: <runs_dir>/<run>__knn_<rep>_k<k>{,_shuff}/ (run-shaped: cv_drug_ranks*.tsv, fold<i>/, compare_real_vs_shuffled/)
#           <runs_dir>/performance_figure/vs_baseline/performance_panels_<metric>.{png,pdf,tsv}
#SBATCH --job-name=baseline_knn
#SBATCH --array=0-2
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 4
#SBATCH --mem=48GB
#SBATCH --time=0-02:00:00
#SBATCH --output=./baseline_knn_%A_%a.out
#SBATCH --error=./baseline_knn_%A_%a.err
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
task=${SLURM_ARRAY_TASK_ID:-0}
if (( task >= ${#runs[@]} )); then echo "array task ${task}: only ${#runs[@]} run(s) given, nothing to do"; exit 0; fi
r="${runs[$task]}"

conda activate metabolic_drug_repurposing
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4} MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK:-4}
status=0
echo "=== ${r}: representations ${reps[*]}, k ${ks[*]}, real + _shuff"
python baseline_signature_knn.py --run-dir "${runs_dir}/${r}" --representation "${reps[@]}" --k "${ks[@]}" || status=1
for rep in "${reps[@]}"; do
    for k in "${ks[@]}"; do
        b="${runs_dir}/${r}__knn_${rep}_k${k}"
        if [[ -s "${b}/cv_drug_ranks.tsv" && -s "${b}_shuff/cv_drug_ranks.tsv" ]]; then
            python compare_runs.py --run real="$b" --run shuffled="${b}_shuff" --control shuffled \
                --out-dir "${b}/compare_real_vs_shuffled" > "${b}/compare.log" 2>&1 || status=1
        else
            echo "WARNING: ${b} or ${b}_shuff incomplete, no comparison"; status=1
        fi
    done
done

# figure: the model rows plus, per signature, each baseline pair that exists (first k only, to keep it readable)
sig() { case "$1" in metchange*) echo mc;; expr*) echo expr;; tf*) echo tf;; esac; }
signame() { case "$1" in mc) echo MetChange;; expr) echo expression;; tf) echo TF;; esac; }
extra=()
for rr in "${runs[@]}"; do
    for rep in "${reps[@]}"; do
        b="${rr}__knn_${rep}_k${ks[0]}"
        [[ -s "${runs_dir}/${b}/cv_drug_ranks.tsv" && -s "${runs_dir}/${b}_shuff/cv_drug_ranks.tsv" ]] || continue
        s=$(sig "$rr")
        what=$([[ "$rep" == sources ]] && echo "kNN, same input" || echo "kNN, full profile")
        extra+=(--extra "${s}=${b}:$(signame "$s") ${what}")
    done
done
if [[ ${#extra[@]} -gt 0 ]]; then
    for m in rank_probc_heldout spec_prob_timedose_heldout rank_prob_heldout; do
        python plot_performance_panels.py --runs-dir "$runs_dir" --out-dir "${runs_dir}/performance_figure/vs_baseline" \
            --metric "$m" "${extra[@]}" || status=1
    done
fi
conda deactivate
exit $status
