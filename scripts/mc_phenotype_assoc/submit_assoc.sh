#!/bin/bash
# Submit the phenotype association for every cell line that HAS usable labels.
# Run this on the login node; it is not a batch script itself.
#
#   bash scripts/mc_phenotype_assoc/submit_assoc.sh                  # all runnable cells
#   cells="HEPG2 A549" bash scripts/mc_phenotype_assoc/submit_assoc.sh
#   dry_run=true bash scripts/mc_phenotype_assoc/submit_assoc.sh     # show, do not submit
#   use_union=true bash scripts/mc_phenotype_assoc/submit_assoc.sh   # SAME phenotype
#                                                                    # table for every cell
#
# USE use_union=true FOR ANY CROSS-CELL-LINE COMPARISON. By default each cell
# line uses its own per-cell phenotype table when one exists and the union table
# otherwise - fine for single-cell runs, but mixed across cell lines, and the two
# tables do not assign identical groups: mapping against the 34,419-compound
# union instead of one cell line's 5,494 lets some external compounds match
# further pert_ids, which changes a few aggregated labels (HEPG2: 144 vs 147
# MMP-down, 475 vs 480 inactive, 27 vs 29 Hallinger-active). Differences of that
# size would otherwise masquerade as cell-line effects.
#
# Why a submit wrapper rather than one batch script per cell:
#   - the set of cell lines with finished limma runs is discovered at run time
#     (check_label_availability.py) instead of hardcoded, so MCF7 or any other
#     line without labels is skipped automatically rather than failing a task;
#   - the Slurm array size has to be known at submit time, so the cell x arm
#     grid is written to a manifest here and each array task reads its own line;
#   - `--output` is expanded by Slurm before the job script runs and cannot
#     contain a shell variable, so passing it here is the only way to put the
#     logs in a fixed absolute directory (scripts/mc_phenotype_assoc/logs).
#
# Cell lines missing their per-cell phenotype tables are reported and skipped:
# run `cell=<CELL> bash prepare_assoc_inputs.sh` for them first.
set -uo pipefail

here="$(cd "$(dirname "$0")" && pwd)"
repo_dir="${repo_dir:-$(cd "$here/../.." && pwd)}"
log_dir="$here/logs"
mkdir -p "$log_dir"

# Resolve the project root through the common config itself, so this logic has
# exactly one implementation (it used to be copied here and into
# check_label_availability.py, and the copies drifted). `cell` is irrelevant for
# the resolution but common insists on it being set.
if [ -z "${labels_base:-}" ]; then
    labels_base="$(cell=_resolve; . "$here/mc_phenotype_assoc_common.sh"; printf '%s' "$labels_base")"
fi
export repo_dir labels_base

py="${py:-python3}"
checker="$here/check_label_availability.py"

echo "repo        : $repo_dir"
echo "labels base : $labels_base"
echo "logs        : $log_dir"
echo
"$py" -I "$checker" --labels-base "$labels_base" ${cells:+--cells "$cells"}
echo

manifest="$log_dir/manifest_$(date +%Y%m%d_%H%M%S).tsv"
: > "$manifest"
n=0
skipped_pheno=""
for arm in mc expr; do
    while read -r c; do
        [ -n "$c" ] || continue
        pheno="$labels_base/results/mc_phenotype_assoc/$c/phenotypes/attene_ramos_mmp.tsv"
        pheno_union="$labels_base/results/mc_phenotype_assoc/phenotypes_ALL/attene_ramos_mmp.tsv"
        if [ ! -f "$pheno" ] && [ ! -f "$pheno_union" ]; then
            case " $skipped_pheno " in *" $c "*) ;; *) skipped_pheno="$skipped_pheno $c";; esac
            continue
        fi
        printf '%s\t%s\n' "$c" "$arm" >> "$manifest"
        n=$((n + 1))
    done <<< "$("$py" -I "$checker" --labels-base "$labels_base" ${cells:+--cells "$cells"} --list "$arm")"
done

if [ -n "$skipped_pheno" ]; then
    echo "SKIPPED (no phenotype tables yet):$skipped_pheno"
    echo "    either per cell:   cell=<CELL> bash $here/prepare_assoc_inputs.sh"
    echo "    or once for all:   cell=ALL   bash $here/prepare_assoc_inputs.sh"
    for c in $skipped_pheno; do
        echo "      missing: $c"
    done
    echo
fi

if [ "$n" -eq 0 ]; then
    echo "nothing to submit" >&2
    rm -f "$manifest"
    exit 1
fi

echo "manifest ($n task(s)): $manifest"
cat "$manifest" | awk '{printf "    [%d] %s %s\n", NR-1, $1, $2}'
echo

export_list="ALL,manifest=$manifest,repo_dir=$repo_dir,labels_base=$labels_base"
if [ "${use_union:-false}" = "true" ]; then
    union_dir="$labels_base/results/mc_phenotype_assoc/phenotypes_ALL"
    if [ ! -f "$union_dir/attene_ramos_mmp.tsv" ]; then
        echo "use_union=true but no union tables at $union_dir" >&2
        echo "build them with: cell=ALL bash $here/prepare_assoc_inputs.sh" >&2
        exit 1
    fi
    export_list="$export_list,pheno_dir=$union_dir"
    echo "using the SAME (union) phenotype table for every cell line: $union_dir"
    echo
fi

cmd=(sbatch --array="0-$((n - 1))"
     --output="$log_dir/mc_pheno_%A_%a.out"
     --export="$export_list"
     "$here/launch_assoc_manifest_slurm.sh")
if [ "${dry_run:-false}" = "true" ]; then
    echo "DRY RUN, would submit:"; printf '    %q ' "${cmd[@]}"; echo
else
    "${cmd[@]}"
fi
