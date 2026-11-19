#!/bin/bash
#SBATCH --job-name=mc_pheno_assoc
#SBATCH --array=0-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/mc_pheno_assoc_%A_%a.out
# Both arms as a 2-task array: 0 = MetChange, 1 = expression baseline.
#   sbatch scripts/mc_phenotype_assoc/launch_assoc_HEPG2_slurm.sh
# Memory: the test holds the selected columns of the label matrix in RAM as
# float32, one matrix per statistic family (~50 MB per matrix for MetChange at
# 5,707 features x ~2,200 perturbations, more for expression at ~20,000 genes).
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
mkdir -p logs

here="scripts/mc_phenotype_assoc"
case "${SLURM_ARRAY_TASK_ID:-0}" in
    0) bash "$here/launch_assoc_metchange_HEPG2.sh" ;;
    1) bash "$here/launch_assoc_expr_HEPG2.sh" ;;
    *) echo "unknown array index" >&2; exit 1 ;;
esac
