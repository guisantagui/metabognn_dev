#!/usr/bin/env bash
# Reproduce every number reported in reports/2026-11-24_prism_gcp_cellpainting.md,
# sections 8 to 15.
#
#   bash scripts/mc_phenotype_assoc/reproduce_report10.sh          # all sections
#   bash scripts/mc_phenotype_assoc/reproduce_report10.sh 12 13    # selected
#
# WHY THIS EXISTS. Those sections were originally run as ad-hoc command lines, so
# the reported numbers were reproducible only by remembering which of
# prism_differential.py's 33 flags had been set. That is nominal reproducibility:
# the script was in git but the RESULT was not, and a changed default would alter
# a published number silently. Section 7 of the same report identified exactly
# this failure for the cross-cell work and it was then repeated a dozen times.
#
# PREREQUISITE, and it is the only one: per-cell scores must already exist under
#   results/mc_phenotype_assoc/<CELL>/prism_viability/per_drug_scores.tsv
#   results/mc_phenotype_assoc/<CELL>/prism_viability_expr/per_drug_scores.tsv
# for A375 A549 HT29 MCF7 PC3. Those come from the cluster:
#   cd scripts/mc_phenotype_assoc && sbatch launch_prism_slurm.sh
# Everything below is local and reads only those ~2 MB files plus the phenotype
# table. Nothing here touches a label matrix.
set -uo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
cd "$repo"
# shellcheck source=python_env.sh
. "$here/python_env.sh"
require_python || exit 1

out="results/mc_phenotype_assoc/cross_cell"
mkdir -p "$out"
PH="results/mc_phenotype_assoc/phenotypes_ALL/prism_viability.tsv"
TGT="data/drug_targets/combined/drug_targets.tsv"
NODES="data/networks/combined/nodes.tsv"
NUC='nucleotide|purine|pyrimidine|folate'
# the hand-made exclusion of 11.3 (used by sections 11-13 as reported)
EXC='^PDE|^TXNRD|^XDH|^GCK$|^HK1$|^EEF2$'
# the principled class list of 14.3/14.5
CLS='^PDE|^PIK3|^PIP[45]K|^SPHK|^DGK|^XDH|^CYP|^UGT|^SULT|^NAT[12]$|^FMO|^EHMT|^EZH|^DNMT|^SETD|^KMT|^PRMT|^COMT$|^NNMT$|^MAO|^ALDH|^ADH[0-9]|^AOX1$|^DBH$|^SIRT|^PARP|^CD38$|^TXNRD|^GSR$|^GPX|^GST|^NQO1$|^SOD[0-9]|^CAT$|^PTGS|^ALOX|^PTGES|^HPGDS$|^TBXAS1$|^NOS[0-9]|^CA[0-9]|^HSD|^SRD5A|^STS$'
D="$py -I $here/prism_differential.py --phenotype $PH"

# Section selection. WANTED holds the script's own arguments; an empty WANTED means
# run everything. `sec` both records the current section and prints its banner only
# when that section is selected - printing unconditionally made a filtered run look
# as though it had executed sections it skipped.
WANTED=("$@")
want() {
    [ ${#WANTED[@]} -eq 0 ] && return 0
    for w in "${WANTED[@]}"; do [ "$w" = "$SEC" ] && return 0; done
    return 1
}
sec() {
    SEC="$1"; shift
    want && printf '\n\n########## report 10, section %s — %s\n' "$SEC" "$*"
    return 0
}

if [ ! -f "results/mc_phenotype_assoc/A549/prism_viability/per_drug_scores.tsv" ]; then
    echo "MISSING per-cell scores; run launch_prism_slurm.sh on the cluster first" >&2
    exit 1
fi

sec 8 "the pooled differential: all drugs, no selection"
want && $D --out "$out/prism_differential.tsv" --top 12

sec 8.2 "the same with scores z-scored within each cell line"
want && $D --standardise-within-cell --out "$out/prism_differential_standardised.tsv" --top 6

sec 9 "MoA-stratified"
want && $D --moa-regex 'antimetabolite|dihydrofolate|thymidylate|ribonucleotide reductase|purine|pyrimidine|nucleoside|nucleotide|folate|inosine monophosphate|dihydroorotate|adenosine deaminase|DNA synthesis|glutamin|amino acid|asparaginase|glycolysis|hexokinase|pyruvate dehydrogenase' \
    --out "$out/prism_differential_moa.tsv" --top 10

sec 10 "target-stratified: any energy+biomass GEM target (90 drugs), then nucleotide only (44)"
if want; then
    $D --targets "$TGT" --nodes "$NODES" --out "$out/prism_differential_targets.tsv" --top 5
    $D --targets "$TGT" --nodes "$NODES" --target-subsystems "$NUC" \
        --out "$out/prism_differential_nucleotide_targets.tsv" --top 0
fi

sec 11 "Jaccard grouping, and the clean 28 after removing the 11.3 contaminants"
if want; then
    for th in 0.5 0.7; do
        $D --targets "$TGT" --nodes "$NODES" --target-subsystems "$NUC" \
            --group-jaccard "$th" --out "$out/prism_diff_jaccard${th}.tsv" --top 4
    done
    for g in "" "--group-jaccard 0.7" "--group-jaccard 0.5"; do
        tag=$(echo "${g:-ungrouped}" | tr -d ' -' | tr '.' '_')
        # shellcheck disable=SC2086
        $D --targets "$TGT" --nodes "$NODES" --target-subsystems "$NUC" \
            --exclude-target-genes "$EXC" --moa-min-drugs 8 $g --top 0 \
            --out "$out/prism_clean28_${tag}.tsv"
    done
    $py -I "$here/prism_differential.py" --phenotype "$PH" --targets "$TGT" --nodes "$NODES" \
        --target-subsystems "$NUC" --group-report 10 --group-report-set dntp_endpoints --top 0 \
        --out "$out/prism_group_report.tsv"
fi

sec 12 "THE HEADLINE: conditional on expression, within the 28 matched drugs"
want && $D --targets "$TGT" --nodes "$NODES" --target-subsystems "$NUC" \
    --exclude-target-genes "$EXC" --expr-run prism_viability_expr --moa-min-drugs 8 --top 0 \
    --out "$out/prism_clean28_conditional.tsv"

sec 13 "specificity: one run per subsystem, drugs and set matched"
if want; then
    declare -a PAIRS=(
      "nucleotide|nucleotide|purine|pyrimidine|folate"
      "purine|^Purine metabolism\$"
      "pyrimidine|^Pyrimidine metabolism\$"
      "folate|^Folate metabolism\$"
      "glycolysis|Glycolysis / Gluconeogenesis"
      "tca|Tricarboxylic acid cycle"
      "oxphos|Oxidative phosphorylation"
      "ppp|Pentose phosphate"
      "pyruvate|^Pyruvate metabolism\$"
      "amino_acid|alanine|arginine|aspartate|cysteine|glutam|glycine|histidine|isoleucine|lysine|methionine|phenylalanine|proline|serine|threonine|tryptophan|tyrosine|valine|urea cycle"
    )
    for spec in "${PAIRS[@]}"; do
        name="${spec%%|*}"; pat="${spec#*|}"
        echo "--- subsystem $name"
        $D --targets "$TGT" --nodes "$NODES" --target-subsystems "$pat" --target-effect inh \
            --exclude-target-genes "$EXC" --expr-run prism_viability_expr \
            --moa-min-drugs 8 --top 0 --out "$out/prism_subsystem_${name}.tsv"
    done
fi

sec 14 "the principled filter stack (14.5 concludes: use the CLASS list alone)"
if want; then
    $py -I "$here/build_target_relevance.py" \
        --out results/mc_phenotype_assoc/metabolite_sets/target_relevance.tsv
    echo "--- class exclusion alone, the 14.5 recommendation"
    $D --targets "$TGT" --nodes "$NODES" --target-subsystems "$NUC" \
        --exclude-target-classes "$CLS" --moa-min-drugs 8 --top 0 \
        --expr-run prism_viability_expr --out "$out/prism_classfilter.tsv"
fi

sec 15 "which cell lines are worth a label run"
want && $py -I "$here/prism_cell_coverage.py" --out "$out/prism_cell_coverage.tsv"

printf '\n\nall outputs in %s\n' "$out"
