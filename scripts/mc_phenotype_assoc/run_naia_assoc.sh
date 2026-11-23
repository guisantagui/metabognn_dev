#!/usr/bin/env bash
# Report 09 section 9 end to end: the Naia SH-SY5Y ATP phenotype against another
# cell line's MetChange labels, with the expression baseline and the controls.
#
# This is NOT a cell-matched test and cannot be made into one: LINCS has 224
# perturbations in SHSY5Y and only 10 of Naia's compounds are among them
# (report 09, 8.1). It fixes the LABEL cell line and varies the PHENOTYPE cell
# line, the converse of the section-3 design - read the direction before reading
# the result.
#
# Usage, from the repository root:
#   bash scripts/mc_phenotype_assoc/run_naia_assoc.sh [CELL]
# CELL is the label cell line, default HEPG2. It needs
# data/LINCS/lincs_unique_perts_<CELL>.tsv and results/labels/<CELL>/.
set -euo pipefail

cell="${1:-HEPG2}"
repo="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo"

sets_dir=results/mc_phenotype_assoc/metabolite_sets
pheno=results/mc_phenotype_assoc/phenotypes_ALL/naia_atp.tsv
mc_out="results/mc_phenotype_assoc/${cell}/naia_atp"
ex_out="results/mc_phenotype_assoc/${cell}/naia_atp_expr"
py() { python3 -I "$@"; }

echo "=== 1/5  phenotype table (groups from the ATP z-score, mapped to the ${cell} panel) ==="
py scripts/mc_phenotype_assoc/build_naia_phenotype.py --cell "$cell" --z-flat 1.0 --out "$pheno"

echo
echo "=== 2/5  are Naia ATP and Attene-Ramos MMP the same readout? (report 09, 9.5) ==="
py scripts/mc_phenotype_assoc/compare_phenotype_tables.py \
    --a results/mc_phenotype_assoc/phenotypes_ALL/attene_ramos_mmp.tsv --a-name attene_MMP \
    --b "$pheno" --b-name naia_ATP --a-case MMP_down --b-case ATP_down || true

echo
echo "=== 3/5  MetChange labels (report 09, 9.1) ==="
py scripts/mc_phenotype_assoc/mc_phenotype_assoc.py \
    --labels "results/labels/${cell}/metchange/rank" \
    --sets "${sets_dir}/metabolite_sets.tsv" --phenotype "$pheno" \
    --mode value --mode sign \
    --case ATP_down --control inactive --extra-contrast ATP_up:inactive \
    --out-dir "$mc_out"

echo
echo "=== 4/5  expression baseline, GPR gene sets (report 09, 9.6) ==="
py scripts/mc_phenotype_assoc/mc_phenotype_assoc.py \
    --labels "results/labels/${cell}/expr" \
    --sets "${sets_dir}/gene_sets.tsv" --phenotype "$pheno" \
    --mode value --mode sign \
    --case ATP_down --control inactive --extra-contrast ATP_up:inactive \
    --out-dir "$ex_out"

echo
echo "=== 5/5  controls: match on MEASURED cytotoxicity, and on the n_sig proxy ==="
echo "--- MetChange ---"
for sv in "glycolysis cyto" "ppp all" "pyruvate cyto" "tca cyto" "energy_currency all"; do
    set -- $sv
    py scripts/mc_phenotype_assoc/covariate_matched_check.py \
        --scores "${mc_out}/per_drug_scores.tsv" --phenotype "$pheno" \
        --set "$1" --variant "$2" --case ATP_down \
        --covariate z_cytotox --covariate n_sig
done
echo "--- expression ---"
for s in oxphos ppp tca glycolysis; do
    py scripts/mc_phenotype_assoc/covariate_matched_check.py \
        --scores "${ex_out}/per_drug_scores.tsv" --phenotype "$pheno" \
        --set "$s" --variant gene --case ATP_down \
        --covariate z_cytotox --covariate n_sig
done
echo
echo "--- robustness of the MetChange headline (report 09, 9.2) ---"
py scripts/mc_phenotype_assoc/covariate_matched_check.py \
    --scores "${mc_out}/per_drug_scores.tsv" --phenotype "$pheno" \
    --set glycolysis --variant cyto --case ATP_down \
    --covariate z_cytotox --covariate n_sig --robustness

echo
echo "done. results in ${mc_out} and ${ex_out}"
