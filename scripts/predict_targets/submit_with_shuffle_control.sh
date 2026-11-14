#!/bin/bash
# Submit a prediction run and its --shuffle-signatures control in parallel, and a comparison job that starts when
# both have ended (afterany: it also runs if one fails, comparing the folds finished in both). For array launchers
# (#SBATCH --array) the dependency on the array job id waits for all its tasks.
# Usage (from this folder, on the login node): bash submit_with_shuffle_control.sh <variant> [<variant> ...]
#        bash submit_with_shuffle_control.sh all   (tf, expr, mc)
#        bash submit_with_shuffle_control.sh <variant> --real-job <jobid>
#   variants: tf expr mc mc_matched mc_matched_auprc mc_matched_noGlobals mc_noGlobals
#             mc_noGlobals_lr1e-3 mc_matched_noGlobals_lr1e-3 mc_noGlobals_lr1e-3_ep100 mc_noGlobals_lr1e-3_seed1
#             tf_noGlobals_lr1e-3 expr_noGlobals_lr1e-3 mc_noGlobals_lr1e-3_min15 tf_noGlobals_lr1e-3_min15
#             expr_noGlobals_lr1e-3_min15 mc_noGlobals_lr1e-3_min15_pat25 mc_noGlobals_lr1e-3_plateau
#             tf_noGlobals_lr1e-3_plateau expr_noGlobals_lr1e-3_plateau (array jobs, one fold per task)
#   (launch_predict_targets_A549_<variant>[_shuff].sh)
#   --real-job: the real run is already queued or running as <jobid> (or finished: give 0); only the control and the
#   comparison are submitted (resubmitting the real launcher would overwrite its out_dir)
# The comparison is written to <real out_dir>/compare_real_vs_shuffled/ and compare_runs_<jobid>.out.
set -euo pipefail
cd "$(dirname "$0")"

results="/home/users/gsantamaria/projects/metabognn_dev/results/predict_targets/A549"
declare -A real_dir=([tf]="${results}/tf_active200" [expr]="${results}/expr" [mc]="${results}/metchange_rank"
                    [mc_matched]="${results}/metchange_rank_matched"
                    [mc_matched_auprc]="${results}/metchange_rank_matched_auprc"
                    [mc_matched_noGlobals]="${results}/metchange_rank_matched_noGlobals"
                    [mc_matched_noGlobals_lr1e-3]="${results}/metchange_rank_matched_noGlobals_lr1e-3"
                    [mc_noGlobals_lr1e-3]="${results}/metchange_rank_noGlobals_lr1e-3"
                    [mc_noGlobals]="${results}/metchange_rank_noGlobals"
                    [mc_noGlobals_lr1e-3_ep100]="${results}/metchange_rank_noGlobals_lr1e-3_ep100"
                    [mc_noGlobals_lr1e-3_seed1]="${results}/metchange_rank_noGlobals_lr1e-3_seed1"
                    [tf_noGlobals_lr1e-3]="${results}/tf_active200_noGlobals_lr1e-3"
                    [expr_noGlobals_lr1e-3]="${results}/expr_noGlobals_lr1e-3"
                    [mc_noGlobals_lr1e-3_min15]="${results}/metchange_rank_noGlobals_lr1e-3_min15"
                    [tf_noGlobals_lr1e-3_min15]="${results}/tf_active200_noGlobals_lr1e-3_min15"
                    [expr_noGlobals_lr1e-3_min15]="${results}/expr_noGlobals_lr1e-3_min15"
                    [mc_noGlobals_lr1e-3_min15_pat25]="${results}/metchange_rank_noGlobals_lr1e-3_min15_pat25"
                    [mc_noGlobals_lr1e-3_plateau]="${results}/metchange_rank_noGlobals_lr1e-3_plateau"
                    [tf_noGlobals_lr1e-3_plateau]="${results}/tf_active200_noGlobals_lr1e-3_plateau"
                    [expr_noGlobals_lr1e-3_plateau]="${results}/expr_noGlobals_lr1e-3_plateau")

out_dir_of() {  # the last out_dir= assignment of a launcher (must match real_dir / real_dir + _shuff)
    grep -E '^out_dir=' "$1" | tail -1
}

submit() {
    local v=$1 real_job=${2:-}  # separate statement: bash expands all of a `local` line before assigning any of it
    local real="launch_predict_targets_A549_${v}.sh" shuff="launch_predict_targets_A549_${v}_shuff.sh"
    echo "== ${v}: $(out_dir_of "$real") | $(out_dir_of "$shuff")"
    local j1 j2 j3
    if [ -n "$real_job" ]; then j1=$real_job; else j1=$(sbatch --parsable "$real"); fi
    j2=$(sbatch --parsable "$shuff")
    local dep="afterany:${j2}"
    [ "$j1" != "0" ] && dep="afterany:${j1}:${j2}"
    j3=$(sbatch --parsable --dependency=${dep} launch_compare_runs.sh \
        --run real="${real_dir[$v]}" --run shuffled="${real_dir[$v]}_shuff" --control shuffled)
    echo "   real ${j1}, shuffled ${j2}, comparison ${j3} (after both)"
}

usage="usage: bash $0 <variant> [<variant> ...] | all | <variant> --real-job <jobid>; variants: ${!real_dir[*]}"
[ $# -eq 0 ] && { echo "$usage"; exit 1; }
if [ "${2:-}" = "--real-job" ]; then
    [ -n "${real_dir[$1]:-}" ] || { echo "$usage"; exit 1; }
    submit "$1" "${3:?--real-job needs a job id (0: already finished)}"
    exit 0
fi
[ "$1" = "all" ] && set -- tf expr mc
for v in "$@"; do  # check all names before submitting anything
    [ -n "${real_dir[$v]:-}" ] || { echo "unknown variant: $v"; echo "$usage"; exit 1; }
done
for v in "$@"; do submit "$v"; done
