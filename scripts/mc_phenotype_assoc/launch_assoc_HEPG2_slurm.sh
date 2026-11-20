#!/bin/bash
#SBATCH --job-name=mc_pheno_HEPG2
#SBATCH --array=0-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=logs/mc_pheno_HEPG2_%A_%a.out
# MetChange / phenotype association for HEPG2.
#
#   cd scripts/mc_phenotype_assoc && sbatch launch_assoc_HEPG2_slurm.sh
#
# SUBMIT FROM scripts/mc_phenotype_assoc: `--output` is expanded by Slurm before
# the job starts and cannot hold a variable, so `logs/` is relative to the submit
# directory. The logs/ directory is kept in git for that reason. For several cell
# lines at once, and to have the log path set absolutely, use submit_assoc.sh
# instead, which also skips cell lines whose labels or phenotype tables are
# missing rather than failing a task.
#
# Can be submitted from anywhere inside the repository: the repo root is located
# by walking up from the submit directory. (Slurm copies the batch script to a
# spool directory before running it, so $0 cannot be relied on to find the repo;
# SLURM_SUBMIT_DIR is where sbatch was invoked.) Export repo_dir to override.
#
# Inputs must exist first, once per cell line:
#   cell=HEPG2 bash scripts/mc_phenotype_assoc/prepare_assoc_inputs.sh
set -uo pipefail

marker="scripts/mc_phenotype_assoc/mc_phenotype_assoc_common.sh"
if [ -z "${repo_dir:-}" ]; then
    for start in "${SLURM_SUBMIT_DIR:-}" "$PWD" "$(dirname "$0")"; do
        [ -n "$start" ] || continue
        d="$(cd "$start" 2>/dev/null && pwd)" || continue
        while [ -n "$d" ] && [ "$d" != "/" ]; do
            if [ -f "$d/$marker" ]; then repo_dir="$d"; break 2; fi
            d="$(dirname "$d")"
        done
    done
fi
if [ -z "${repo_dir:-}" ]; then
    echo "cannot locate the repository root (looked for $marker upward from" >&2
    echo "  SLURM_SUBMIT_DIR=${SLURM_SUBMIT_DIR:-unset}, PWD=$PWD)." >&2
    echo "Export repo_dir=/path/to/metabognn_dev and resubmit." >&2
    exit 1
fi

export cell=HEPG2 repo_dir
exec bash "$repo_dir/scripts/mc_phenotype_assoc/launch_assoc_slurm_generic.sh"
