#!/bin/bash
# MetChange labels vs the measured mitochondrial phenotypes, all three modalities.
# CELL-AGNOSTIC: the cell line comes from $cell, which the slurm launchers export
# (launch_assoc_HEPG2_slurm.sh / launch_assoc_A549_slurm.sh). Run directly with
#   cell=HEPG2 bash launch_assoc_metchange.sh
#
#   value  magnitude  (z, auroc)
#   sign   direction  (sign_mean, sign_dev)
#   sig    significance (sig_dev, sig_signed_dev; needs qvals.csv)
#
# Attene-Ramos (membrane potential, 147 MMP-down / 480 inactive / 75 MMP-up /
# 24 cytotoxic) and Hallinger (OCR, 29 active / 86 inactive, and 19 ETC
# inhibitors), plus the responsiveness-stratified check for each.
#
#   cell=HEPG2 bash launch_assoc_metchange.sh                       # everything
#   cell=A549 modes="--mode sign" bash launch_assoc_metchange.sh    # one modality
#
# Labels, phenotypes, sets and outputs all live under $scratch_proj_dir, which
# mc_phenotype_assoc_common.sh points at scratch on the cluster and at the repo
# on the laptop. Phenotype tables are per cell line; metabolite sets are shared.
set -uo pipefail
. "$(dirname "$0")/mc_phenotype_assoc_common.sh"
provenance
require_inputs || exit 1

# value and sign are defined for every perturbation; sig is not - with no
# significant feature sig_dev is identically 0, so pooling those drugs makes the
# test detect responsiveness rather than pathway enrichment (report 07, 11.6).
# Hence two passes: value+sign on everything, sig on the responsive subset.
modes="${modes:---mode value --mode sign}"
sig_modes="${sig_modes:---mode sig --min-nsig 1}"

run_assoc attene_allmodes "$mc_labels" "$met_sets" \
    "$pheno_dir/attene_ramos_mmp.tsv" MMP_down inactive \
    $modes --extra-contrast MMP_up:inactive --extra-contrast cytotox:inactive
run_strat attene_allmodes MMP_down
run_strat attene_allmodes MMP_up
run_matched attene_allmodes MMP_down z
run_matched attene_allmodes MMP_down sign_dev
run_matched attene_allmodes MMP_up z

run_assoc hallinger_allmodes "$mc_labels" "$met_sets" \
    "$pheno_dir/hallinger_ocr.tsv" active inactive $modes
run_strat hallinger_allmodes active

if [ -f "$pheno_dir/hallinger_ocr_bymech.tsv" ]; then
    run_assoc hallinger_etci "$mc_labels" "$met_sets" \
        "$pheno_dir/hallinger_ocr_bymech.tsv" ETCi inactive \
        $modes --extra-contrast uncoupler:inactive
fi

# significance modality, responsive perturbations only
run_assoc attene_sig "$mc_labels" "$met_sets" \
    "$pheno_dir/attene_ramos_mmp.tsv" MMP_down inactive \
    $sig_modes --extra-contrast MMP_up:inactive
run_assoc hallinger_sig "$mc_labels" "$met_sets" \
    "$pheno_dir/hallinger_ocr.tsv" active inactive $sig_modes

# value/sign on the MetChange-responsive subset, the sensitivity analysis of 11.4
run_assoc attene_responsive "$mc_labels" "$met_sets" \
    "$pheno_dir/attene_ramos_mmp.tsv" MMP_down inactive \
    $modes --min-nsig 1

# declared negative controls (report 08, 3.1): glycan / GAG subsystems chosen on
# biology alone. Four of five were null on HEPG2 and one (O-glycan) was not.
if [ -f "$negctrl_sets" ]; then
    run_assoc attene_negctrl "$mc_labels" "$negctrl_sets" \
        "$pheno_dir/attene_ramos_mmp.tsv" MMP_down inactive \
        --mode value --mode sign --extra-contrast MMP_up:inactive
    run_matched attene_negctrl MMP_down z
fi

say "done; results under $out_base"
