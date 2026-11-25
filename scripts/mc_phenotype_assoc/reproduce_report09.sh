#!/usr/bin/env bash
# Reproduce the numbers in reports/2026-11-23_cross_cell_specificity.md (report 09).
#
#   bash scripts/mc_phenotype_assoc/reproduce_report09.sh          # sections 3,4,5
#   bash scripts/mc_phenotype_assoc/reproduce_report09.sh 4        # one section
#
# WHY. The HepG2 MMP association itself was always scripted - launch_assoc_HEPG2_slurm.sh
# runs the MetChange and expression arms with their stratified and responsiveness-matched
# checks. What had no launcher was the CROSS-CELL comparison of section 3 and the
# MetChange-versus-expression head-to-head of section 4, both of which apply the
# common-drug restriction that the per-cell launchers do not. Section 4 contains the
# project's single most load-bearing number for the central thesis - MetChange 0.365
# against expression 0.500 under the same matched design - and it was reproducible only
# from an ad-hoc command line.
#
# PREREQUISITE: per-cell association outputs for the eight cells with MetChange labels,
# i.e. results/mc_phenotype_assoc/<CELL>/{attene_allmodes,attene_expr,attene_negctrl}/.
# Those come from the cluster:
#   cd scripts/mc_phenotype_assoc && bash submit_assoc.sh
# Everything here is local and reads only per_drug_scores.tsv files.
#
# For report 09 section 9 / 9.6 (the Naia ATP arm, including the MetChange-versus-
# expression comparison there) use run_naia_assoc.sh instead - that one always existed.
set -uo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
cd "$repo"
# shellcheck source=python_env.sh
. "$here/python_env.sh"
require_python || exit 1

out="results/mc_phenotype_assoc/cross_cell"
mkdir -p "$out"
CC="$py -I $here/cross_cell_compare.py"
CELLS="HEPG2 A549 A375 HA1E HCC515 HT29 PC3 VCAP"

WANTED=("$@")
want() { [ ${#WANTED[@]} -eq 0 ] && return 0; for w in "${WANTED[@]}"; do [ "$w" = "$SEC" ] && return 0; done; return 1; }
sec() { SEC="$1"; shift; want && printf '\n\n########## report 09, section %s — %s\n' "$SEC" "$*"; return 0; }

if [ ! -f "results/mc_phenotype_assoc/HEPG2/attene_allmodes/per_drug_scores.tsv" ]; then
    echo "MISSING per-cell scores; run submit_assoc.sh on the cluster first" >&2
    exit 1
fi

sec 3 "cross-cell, MetChange, MMP-down: HEPG2 only"
if want; then
    for sv in "tca mito_nohub" "tca_endpoints cyto_nohub" "tca cyto" "oxphos mito" "energy_currency mito"; do
        set -- $sv
        $CC --cells "$CELLS" --run attene_allmodes --set "$1" --variant "$2" --case MMP_down \
            --out "$out/r09s3_${1}_${2}_MMP_down.tsv"
    done
fi

sec 3.1 "cross-cell, MMP-up: HCC515 survives too, so the specificity is MMP-down only"
if want; then
    $CC --cells "$CELLS" --run attene_allmodes --set energy_currency --variant mito --case MMP_up \
        --out "$out/r09s3.1_energy_currency_mito_MMP_up.tsv"
    $CC --cells "$CELLS" --run attene_allmodes --set tca --variant mito_nohub --case MMP_up \
        --out "$out/r09s3.1_tca_mito_nohub_MMP_up.tsv"
fi

sec 4 "THE HEAD-TO-HEAD: expression collapses to 0.500 in HEPG2 where MetChange holds at 0.365"
if want; then
    echo "--- MetChange arm"
    $CC --cells "$CELLS" --run attene_allmodes --set tca --variant mito_nohub --case MMP_down \
        --out "$out/r09s4_metchange_tca_mito_nohub.tsv"
    echo "--- expression arm, same drugs and same matched design"
    $CC --cells "$CELLS" --run attene_expr --set tca --variant gene --case MMP_down \
        --out "$out/r09s4_expr_tca_gene.tsv"
    $CC --cells "$CELLS" --run attene_expr --set pyruvate --variant gene --case MMP_up \
        --out "$out/r09s4_expr_pyruvate_gene_MMP_up.tsv"
fi

sec 5 "the O-glycan negative control is not cell-specific either"
want && $CC --cells "$CELLS" --run attene_negctrl --set negctrl_o_glycan --variant all --case MMP_down \
    --out "$out/r09s5_negctrl_o_glycan.tsv"

printf '\n\nall outputs in %s\n' "$out"
printf 'for report 09 sections 9 and 9.6 (Naia ATP) use: bash %s/run_naia_assoc.sh\n' "$here"
