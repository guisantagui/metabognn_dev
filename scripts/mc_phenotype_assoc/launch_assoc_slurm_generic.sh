#!/bin/bash
#SBATCH --job-name=mc_pheno_assoc
#SBATCH --array=0-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/mc_pheno_assoc_%x_%A_%a.out
# Generic body shared by the per-cell slurm launchers. Do not sbatch this
# directly: it takes the cell line from $cell, which the per-cell wrappers set
# and export (launch_assoc_HEPG2_slurm.sh, launch_assoc_A549_slurm.sh).
#
# Array: task 0 = MetChange arm, task 1 = expression baseline.
#
# Memory: the test holds the selected columns of the label matrix in RAM as
# float32, one matrix per statistic family (~50 MB for MetChange at 5,707
# features x ~2,200 perturbations, ~4x that for expression at ~20,000 genes).
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
mkdir -p logs

: "${cell:?this script must be called by a per-cell wrapper that exports cell}"
export cell
here="scripts/mc_phenotype_assoc"

case "${SLURM_ARRAY_TASK_ID:-0}" in
    0) bash "$here/launch_assoc_metchange.sh" ;;
    1) bash "$here/launch_assoc_expr.sh" ;;
    *) echo "unknown array index ${SLURM_ARRAY_TASK_ID:-}" >&2; exit 1 ;;
esac
