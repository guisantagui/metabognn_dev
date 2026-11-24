#!/bin/bash
# Score the PRISM drug set against ONE cell line's MetChange labels.
#
# Called per array task by launch_prism_slurm.sh; `cell` must be exported.
# Runs standalone too:
#     cell=A549 bash scripts/mc_phenotype_assoc/launch_prism_cell.sh
#
# WHAT THIS PRODUCES AND WHAT TO IGNORE. The output that matters is
# `per_drug_scores.tsv` - one set score per drug per set x variant x statistic -
# which is the predictor side of the between-cell differential. It is about 3 MB,
# so **this is the only file that needs downloading**; the GB-scale label matrices
# stay on the cluster.
#
# `set_results.tsv` contains the sensitive-vs-resistant group contrast. **That is
# the confounded naive test and must not be read as a result** (report 10, 1.1):
# the phenotype is viability, the confounder of reports 07-09 is viability, and you
# cannot match on cytotoxicity when cytotoxicity is the phenotype. The real test is
# `prism_differential.py`, which differences the confounder out.
#
# Because the group contrast is not the result, its resampling is wasted work and
# is turned down to the minimum here (the differential computes its own null by
# permuting drug identity). Raise n_perm / n_boot / n_random_sets only if you
# actually want the per-cell contrast for its own sake.
set -uo pipefail

: "${cell:?set cell, e.g. cell=A549}"
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=mc_phenotype_assoc_common.sh
. "$here/mc_phenotype_assoc_common.sh"

pheno="${prism_pheno:-$assoc_base/phenotypes_ALL/prism_viability.tsv}"
n_perm="${n_perm:-200}"
n_boot="${n_boot:-200}"
n_random_sets="${n_random_sets:-10}"

provenance
require_python || exit 1
printf 'PRISM table      : %s\n' "$pheno"
printf 'NOTE             : read per_drug_scores.tsv; set_results.tsv is the confounded naive test\n\n'

if [ ! -f "$pheno" ]; then
    cat <<MSG >&2
MISSING PRISM phenotype table: $pheno

Build it on the machine that holds the PRISM download (it is a few hundred kB)
and copy that one file across - the PRISM CSVs themselves do not need to be here:

    python3 -I scripts/mc_phenotype_assoc/build_prism_phenotype.py \\
        --out results/mc_phenotype_assoc/phenotypes_ALL/prism_viability.tsv
    scp results/mc_phenotype_assoc/phenotypes_ALL/prism_viability.tsv \\
        <host>:$assoc_base/phenotypes_ALL/
MSG
    exit 1
fi
if [ ! -d "$mc_labels" ]; then
    echo "MISSING MetChange labels for $cell: $mc_labels" >&2
    exit 1
fi

run_assoc prism_viability "$mc_labels" "$met_sets" "$pheno" sensitive resistant \
    --mode value --mode sign

out="$out_base/prism_viability/per_drug_scores.tsv"
if [ -f "$out" ]; then
    printf '\n=== done: %s (%s)\n' "$out" "$(du -h "$out" | cut -f1)"
    printf 'Download just this file for each cell line, then run locally or here:\n'
    printf '  python3 -I scripts/mc_phenotype_assoc/prism_differential.py \\\n'
    printf '      --phenotype %s \\\n' "$pheno"
    printf '      --scores-root %s\n' "$assoc_base"
else
    echo "FAILED: no $out" >&2
    exit 1
fi
