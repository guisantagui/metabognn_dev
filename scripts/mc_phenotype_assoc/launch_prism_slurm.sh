#!/bin/bash
#SBATCH --job-name=prism_mc
#SBATCH --array=0-4
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/prism_mc_%A_%a.out
# PRISM viability x MetChange, one array task per PRISM-covered LINCS core line.
#
#   cd scripts/mc_phenotype_assoc && sbatch launch_prism_slurm.sh
#
# SUBMIT FROM scripts/mc_phenotype_assoc: `--output` is expanded by Slurm before
# the job starts and cannot hold a variable, so `logs/` is relative to the submit
# directory (which is why logs/ is kept in git).
#
# The five lines are A375, A549, HT29, MCF7, PC3 - the LINCS cores present in the
# PRISM secondary screen. **HEPG2, VCAP and HCC515 are absent from PRISM** and
# HA1E is absent by construction, so the HepG2 results of reports 07-09 cannot be
# extended here, only tested elsewhere (report 10, section 1).
#
# Each task needs that cell line's MetChange labels under
# <scratch_proj_dir>/results/labels/<cell>/metchange/rank, plus the single small
# phenotype table at <scratch_proj_dir>/results/mc_phenotype_assoc/phenotypes_ALL/
# prism_viability.tsv (build it where the PRISM download is and scp it over).
#
# A missing cell line fails only its own task, so run what you have and re-run the
# rest later; prism_differential.py uses whichever cells it finds and says which.
#
# Each task runs BOTH arms: MetChange labels with the metabolite sets, and the
# expression baseline with the GPR gene sets. The expression arm is what lets
# prism_differential.py ask whether MetChange adds information the transcriptome
# lacks - the project's actual claim - and these are the first five non-HepG2 cell
# lines where that can be asked (report 10, 10.5). For one arm only:
#   arm=metchange sbatch launch_prism_slurm.sh      # or arm=expr
#
# Override the cell list with: cells="A549 PC3" sbatch --array=0-1 launch_prism_slurm.sh
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

cells="${cells:-A375 A549 HT29 MCF7 PC3}"
read -r -a _cells <<< "$cells"
idx="${SLURM_ARRAY_TASK_ID:-0}"
if [ "$idx" -ge "${#_cells[@]}" ]; then
    echo "array index $idx beyond the cell list ($cells)" >&2
    exit 1
fi

export cell="${_cells[$idx]}" repo_dir arm="${arm:-both}"
echo "array task $idx -> cell=$cell"
exec bash "$repo_dir/scripts/mc_phenotype_assoc/launch_prism_cell.sh"
