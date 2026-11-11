#!/bin/bash
# Submit a prediction run and its --shuffle-signatures control in parallel, and a comparison job that starts when
# both have ended (afterany: it also runs if one fails, comparing the folds finished in both).
# Usage (from this folder, on the login node): bash submit_with_shuffle_control.sh <tf|expr|mc|all>
#        bash submit_with_shuffle_control.sh <tf|expr|mc> --real-job <jobid>
#   --real-job: the real run is already queued or running as <jobid> (or finished: give 0); only the control and the
#   comparison are submitted (resubmitting the real launcher would overwrite its out_dir)
# The comparison is written to <real out_dir>/compare_real_vs_shuffled/ and compare_runs_<jobid>.out.
set -euo pipefail
cd "$(dirname "$0")"

results="/home/users/gsantamaria/projects/metabognn_dev/results/predict_targets/A549"
declare -A real_dir=([tf]="${results}/tf_active200" [expr]="${results}/expr" [mc]="${results}/metchange_rank")

out_dir_of() {  # the last out_dir= assignment of a launcher (must match real_dir / real_dir + _shuff)
    grep -E '^out_dir=' "$1" | tail -1
}

submit() {
    local v=$1 real_job=${2:-} real="launch_predict_targets_A549_${v}.sh" shuff="launch_predict_targets_A549_${v}_shuff.sh"
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

case "${1:-}" in
    tf|expr|mc)
        if [ "${2:-}" = "--real-job" ]; then submit "$1" "${3:?--real-job needs a job id (0: already finished)}"
        else submit "$1"; fi ;;
    all) for v in tf expr mc; do submit "$v"; done ;;
    *) echo "usage: bash $0 <tf|expr|mc|all>"; exit 1 ;;
esac
