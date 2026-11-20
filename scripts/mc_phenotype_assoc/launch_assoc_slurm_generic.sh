#!/bin/bash
# Generic body shared by the per-cell slurm launchers. Not for direct sbatch:
# the cell line and the repo root come from the wrapper, which exports them
# (launch_assoc_HEPG2_slurm.sh, launch_assoc_A549_slurm.sh).
#
# Array: task 0 = MetChange arm, task 1 = expression baseline.
#
# Memory: the test holds the selected columns of the label matrix in RAM as
# float32, one matrix per statistic family (~50 MB for MetChange at 5,707
# features x ~2,200 perturbations, ~4x that for expression at ~20,000 genes).
set -uo pipefail

: "${cell:?call this through a per-cell wrapper that exports cell}"
: "${repo_dir:?call this through a per-cell wrapper that exports repo_dir}"
export cell repo_dir

here="$repo_dir/scripts/mc_phenotype_assoc"
cd "$repo_dir" || exit 1

case "${SLURM_ARRAY_TASK_ID:-0}" in
    0) bash "$here/launch_assoc_metchange.sh" ;;
    1) bash "$here/launch_assoc_expr.sh" ;;
    *) echo "unknown array index ${SLURM_ARRAY_TASK_ID:-}" >&2; exit 1 ;;
esac
