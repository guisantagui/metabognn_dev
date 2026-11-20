#!/bin/bash
#SBATCH --job-name=mc_pheno
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/mc_pheno_%A_%a.out
# Array body for the multi-cell submission. Do not sbatch directly: submit_assoc.sh
# builds the manifest, sizes the array and passes an absolute --output (which
# overrides the directive above, since Slurm expands --output before the job
# runs and it cannot contain a variable).
#
# Each array task handles one (cell line, arm) pair read from its manifest line.
set -uo pipefail

: "${manifest:?submit through submit_assoc.sh}"
: "${repo_dir:?submit through submit_assoc.sh}"
idx="${SLURM_ARRAY_TASK_ID:-0}"

line="$(awk -v i="$((idx + 1))" 'NR==i' "$manifest")"
if [ -z "$line" ]; then
    echo "no manifest line $idx in $manifest" >&2
    exit 1
fi
cell="$(printf '%s' "$line" | cut -f1)"
arm="$(printf '%s' "$line" | cut -f2)"
export cell repo_dir

echo "array task $idx -> cell=$cell arm=$arm   (manifest $manifest)"
cd "$repo_dir" || exit 1

case "$arm" in
    mc)   exec bash "$repo_dir/scripts/mc_phenotype_assoc/launch_assoc_metchange.sh" ;;
    expr) exec bash "$repo_dir/scripts/mc_phenotype_assoc/launch_assoc_expr.sh" ;;
    *)    echo "unknown arm '$arm'" >&2; exit 1 ;;
esac
