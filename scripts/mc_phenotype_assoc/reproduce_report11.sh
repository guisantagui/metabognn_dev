#!/usr/bin/env bash
# Reproduce reports/2026-11-25_nci60_exchange_flux_feasibility.md, sections 8 to 10.
#
#   bash scripts/mc_phenotype_assoc/reproduce_report11.sh          # all sections
#   bash scripts/mc_phenotype_assoc/reproduce_report11.sh 9 10     # selected
#   cell=A549 bash scripts/mc_phenotype_assoc/reproduce_report11.sh
#
# WHAT IT ANSWERS
#   8   How many independent dimensions does a label matrix have, and is that
#       unusual? (compared against the expression arm of the same cells)
#   9   The producer-average null: can MetChange beat the mean expression of the
#       enzymes that produce each metabolite - same genes, no LP, no stoichiometry?
#   10  Does MetChange's advantage survive removing the dominant global component?
#   11  The four-way table, all arms on identical (contrast, set, variant, statistic)
#       keys. Needs 9 and 10 to have run. Read `rand`, not the parametric p.
#
# PREREQUISITE: local label matrices for the cell line, BOTH arms:
#   results/labels/<CELL>/metchange/rank/log2FC.csv
#   results/labels/<CELL>/expr/log2FC.csv
# As of 2026-11-25 that is HEPG2 and A549 only. HEPG2 is the default because it is
# the one cell line where MetChange demonstrably beat expression (report 09, 4),
# which is the result sections 9 and 10 are trying to break.
#
# RUNTIME: ~15 min for all three sections on a laptop. Section 9 loads the GEM
# (7 s) and both label matrices; nothing here needs the cluster.
set -uo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
cd "$repo"
# shellcheck source=python_env.sh
. "$here/python_env.sh"
require_python || exit 1

cell="${cell:-HEPG2}"

# WHERE THE LABELS LIVE. On the cluster they are on scratch, on a laptop in the repo.
# mc_phenotype_assoc_common.sh already resolves this as `labels_base`; ask it rather than
# assuming repo-relative paths, which is what every other launcher here does and what an
# earlier version of this script failed to do (it would have reported "needs BOTH label
# arms locally" on the cluster while the matrices sat on scratch).
if [ -z "${labels_base:-}" ]; then
    labels_base="$(cell="$cell"; . "$here/mc_phenotype_assoc_common.sh" >/dev/null 2>&1; printf '%s' "$labels_base")"
fi
: "${labels_base:=$repo}"
export LABELS_BASE="$labels_base"      # read by feature_dimensionality.py / build_derived_labels.py
echo "repo        : $repo"
echo "labels base : $labels_base"
# The ASSAY cell line of each phenotype. The phenotype never changes cell line when
# `cell` changes - only the LABELS do. So running a HepG2-measured phenotype against
# another cell's labels is the cross-cell CONTROL of report 08 section 8.3, not a
# replication: a null there means the association is that cell's metabolism rather
# than a compound property. Printed below because this reading has been inverted
# once already (report 11, 12.0).
pheno_cell="HEPG2"
sets_dir="$labels_base/results/mc_phenotype_assoc/metabolite_sets"
PH="$labels_base/results/mc_phenotype_assoc/phenotypes_ALL/attene_ramos_mmp.tsv"
out_root="$labels_base/results/mc_phenotype_assoc/${cell}"
mkdir -p "$out_root"

if [ ! -f "$labels_base/results/labels/${cell}/metchange/rank/log2FC.csv" ] ||
   [ ! -f "$labels_base/results/labels/${cell}/expr/log2FC.csv" ]; then
    echo "ERROR: ${cell} needs BOTH label arms locally:" >&2
    echo "  $labels_base/results/labels/${cell}/metchange/rank/log2FC.csv" >&2
    echo "  $labels_base/results/labels/${cell}/expr/log2FC.csv" >&2
    echo "available: $(ls "$labels_base/results/labels" 2>/dev/null | tr '\n' ' ')" >&2
    exit 1
fi

# Are the two arms co-extensive? If not, every arm must be cut down to the shared
# perturbations or the comparison is between different drug sets, not between
# representations. HEPG2 matches exactly (15,032 both); A549 does not (44,396 vs
# 7,636, sharing 7,633). Decided here, once, rather than left to the caller.
n_mc=$(head -1 "$labels_base/results/labels/${cell}/metchange/rank/log2FC.csv" | tr ',' '\n' | wc -l | tr -d ' ')
n_ex=$(head -1 "$labels_base/results/labels/${cell}/expr/log2FC.csv" | tr ',' '\n' | wc -l | tr -d ' ')
if [ "$n_mc" != "$n_ex" ]; then
    RESTRICT="--restrict-to-shared"
    NO_NSIG="--no-nsig"
    MC_ARM="$labels_base/results/labels/${cell}/metchange_shared"
    echo "NOTE: ${cell} arms are NOT co-extensive (metchange $((n_mc-1)) vs expr $((n_ex-1)))."
    echo "      Restricting every arm to the shared perturbations."
else
    RESTRICT=""
    NO_NSIG=""
    MC_ARM="$labels_base/results/labels/${cell}/metchange/rank"
    echo "NOTE: ${cell} arms are co-extensive ($((n_mc-1)) perturbations); no restriction needed."
fi
echo

if [ "$cell" = "$pheno_cell" ]; then
    echo "DESIGN: phenotype measured in ${pheno_cell}, labels from ${cell} -> MATCHED."
    echo "        A positive association is the result under test."
else
    echo "DESIGN: phenotype measured in ${pheno_cell}, labels from ${cell} -> CROSS-CELL CONTROL."
    echo "        This is report 08 section 8.3's control, NOT a replication of the ${pheno_cell} run."
    echo "        A NULL here is the control PASSING: it says the association is ${pheno_cell}"
    echo "        metabolism rather than a compound property. A positive result would be bad news."
fi
echo

WANTED=("$@")
want() {
    [ ${#WANTED[@]} -eq 0 ] && return 0
    for w in "${WANTED[@]}"; do [ "$w" = "$SEC" ] && return 0; done
    return 1
}
sec() {
    SEC="$1"; shift
    want || return 1
    echo; echo "############ section ${SEC}: $* ############"; echo
    return 0
}

# One association run. $1 = label dir, $2 = sets file, $3 = output subdir.
# NO_NSIG is set when the arms had to be restricted: the carried-over q-values and
# n_sig diagnostics index the unrestricted perturbation set, so they are absent from
# the derived dirs by design and must not be looked for. See build_derived_labels.py.
# --mode sign is NOT optional here. `--mode value` gives `z` and `auroc`, both of which
# are computed WITHIN a perturbation (z subtracts that perturbation's mean over all
# features; auroc is a within-column ECDF), so both are RELATIVE statistics and neither
# carries an absolute direction. `--mode sign` adds `sign_mean`, whose baseline is 'raw'
# and which is therefore the only uncentred directional statistic in the pipeline:
# positive = the set's members moved UP (harder to produce), negative = down.
#
# Running value alone is how the TCA/MMP direction was misreported twice (report 11
# 10.3b): first as "easier" by reading z as a log2FC, then as "harder" from a raw mean of
# a heavy-tailed distribution. `sign_mean` settles it in one number - and for
# `tca all_nohub` it is AUROC 0.523 at random-set p 0.995, i.e. no absolute effect at all.
assoc() {
    "$py" -I "$here/mc_phenotype_assoc.py" \
        --labels "$1" --sets "$2" --phenotype "$PH" \
        --mode value --mode sign ${NO_NSIG:-} \
        --case MMP_down --control inactive --extra-contrast MMP_up:inactive \
        --out-dir "$3"
}

# The report 09 section 4 matched comparison, on whichever scores it is given.
matched() {
    for sv in "$@"; do
        set -- $sv
        "$py" -I "$here/covariate_matched_check.py" \
            --scores "${SCORES}/per_drug_scores.tsv" --phenotype "$PH" \
            --set "$1" --variant "$2" --case MMP_down \
            --covariate z_cytotox --covariate n_sig || true
    done
}

############################################################################
if sec 8 "effective dimensionality of every representation"; then
    "$py" -I "$here/feature_dimensionality.py" --cells "$cell"
fi

############################################################################
if sec 9 "the producer-average null (report 11, section 6)"; then
    echo "--- 9.1 build the producer-average matrix ---"
    "$py" -I "$here/build_derived_labels.py" --cell "$cell" --mode producer_average $RESTRICT

    echo
    echo "--- 9.2 its dimensionality, against MetChange and raw expression ---"
    echo "    (if the producer average is NOT collapsed, the LP is what collapses it)"
    "$py" -I - <<'PY'
import sys, os
sys.path.insert(0, os.path.join(os.getcwd(), 'scripts/mc_phenotype_assoc'))
from feature_dimensionality import load_matrix, spectrum
cell = os.environ.get('cell', 'HEPG2')
for lab, p in [('MetChange (LP)', os.path.join(os.environ.get('LABELS_BASE','.'), 'results/labels/%s/metchange/rank/log2FC.csv' % cell)),
               ('producer-average (no LP)', os.path.join(os.environ.get('LABELS_BASE','.'), 'results/labels/%s/producer_average/log2FC.csv' % cell)),
               ('expression (genes)', os.path.join(os.environ.get('LABELS_BASE','.'), 'results/labels/%s/expr/log2FC.csv' % cell))]:
    if not os.path.exists(p):
        continue
    X, _, _ = load_matrix(p, max_rows=8000)
    s = spectrum(X)
    print('  %-26s %6dx%-6d  PC1 %5.1f%%  PR %7.2f  eig@90%% %5d  median r %7.4f'
          % (lab, s['n_features'], s['n_perturbations'], s['pc1'],
             s['participation_ratio'], s['n_eig_90'], s['median_r']))
PY

    echo
    echo "--- 9.3 association, producer-average arm ---"
    assoc "results/labels/${cell}/producer_average" \
          "${sets_dir}/metabolite_sets.tsv" "${out_root}/attene_prodavg"
    echo
    echo "--- 9.4 the same matched comparison report 09 section 4 ran on MetChange ---"
    SCORES="${out_root}/attene_prodavg"
    matched "tca cyto" "oxphos mito" "glycolysis cyto" "dntp_endpoints all"
fi

############################################################################
if sec 10 "does the advantage survive removing the global component?"; then
    echo "--- 10.1 build the residual matrix ---"
    "$py" -I "$here/build_derived_labels.py" --cell "$cell" --mode global_residual $RESTRICT

    echo
    echo "--- 10.2 association, global-residual arm ---"
    assoc "results/labels/${cell}/global_residual" \
          "${sets_dir}/metabolite_sets.tsv" "${out_root}/attene_resid"
    echo
    echo "--- 10.3 matched comparison, residual arm ---"
    SCORES="${out_root}/attene_resid"
    matched "tca cyto" "oxphos mito" "glycolysis cyto" "dntp_endpoints all"

    echo
    echo "--- 10.4 the two reference arms, for the four-way table ---"
    for arm in metchange expr; do
        if [ "$arm" = metchange ]; then
            if [ -n "$RESTRICT" ] && [ ! -f "$MC_ARM/log2FC.csv" ]; then
                "$py" -I "$here/build_derived_labels.py" --cell "$cell" --mode metchange_shared
            fi
            ldir="$MC_ARM"; sf="${sets_dir}/metabolite_sets.tsv"
            odir="${out_root}/attene_mc_ref"
        else
            ldir="$labels_base/results/labels/${cell}/expr"; sf="${sets_dir}/gene_sets.tsv"
            odir="${out_root}/attene_expr_ref"
            if [ -n "$RESTRICT" ]; then
                ex_arm="$labels_base/results/labels/${cell}/expr_shared"
                if [ ! -f "$ex_arm/log2FC.csv" ]; then
                    "$py" -I - <<PYEOF
import sys, os
sys.path.insert(0, "$here")
from build_derived_labels import shared_columns, restrict_matrix, write_labels, label_path
cols = shared_columns("$cell")
M, ids, kept = restrict_matrix("$cell", "expr", cols)
write_labels(M, ids, kept, "$ex_arm", os.path.dirname(label_path("$cell", "expr")))
PYEOF
                fi
                ldir="$ex_arm"
            fi
        fi
        echo "  [$arm]"
        assoc "$ldir" "$sf" "$odir"
    done
fi

############################################################################
if sec 11 "the four-way table"; then
    echo "--- z: the set RELATIVE to the rest of the matrix (what reports 09-11 quote) ---"
    "$py" -I "$here/compare_arms.py" --base "$out_root" --statistic z
    echo
    echo "--- auroc: the same, rank-based, so robust to the heavy tails ---"
    "$py" -I "$here/compare_arms.py" --base "$out_root" --statistic auroc
    echo
    echo "--- sign_mean: UNCENTRED, the only absolute DIRECTION available ---"
    echo "    median_case > 0 means the set moved UP (harder to produce) in the cases."
    "$py" -I "$here/compare_arms.py" --base "$out_root" --statistic sign_mean
    echo
    echo "--- sign_dev: direction relative to the perturbation's own mean sign ---"
    "$py" -I "$here/compare_arms.py" --base "$out_root" --statistic sign_dev
fi

echo
echo "============================================================"
echo "done. four-way comparison, same phenotype and design:"
echo "  MetChange        ${out_root}/attene_mc_ref/set_results.tsv"
echo "  producer-average ${out_root}/attene_prodavg/set_results.tsv"
echo "  global residual  ${out_root}/attene_resid/set_results.tsv"
echo "  expression       ${out_root}/attene_expr_ref/set_results.tsv"
echo
echo "Read the 'rand' column, not the parametric p: with features correlated at"
echo "r=0.97 the parametric null is not the relevant one (report 11, 8.1)."
