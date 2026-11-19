#!/bin/bash
# MetChange labels vs the measured mitochondrial phenotypes, all three modalities.
#
#   value  magnitude  (z, auroc)
#   sign   direction  (sign_mean, sign_dev)
#   sig    significance (sig_dev, sig_signed_dev; needs qvals.csv)
#
# Attene-Ramos (membrane potential, 147 MMP-down / 480 inactive / 75 MMP-up /
# 24 cytotoxic) and Hallinger (OCR, 29 active / 86 inactive, and 19 ETC
# inhibitors), plus the responsiveness-stratified check for each.
#
#   bash launch_assoc_metchange_HEPG2.sh            # everything
#   modes="--mode sig" bash launch_assoc_metchange_HEPG2.sh   # one modality
#
# On the cluster the labels are read from
# /scratch/users/gsantamaria/test_large_files/metabognn_dev/results/labels/<cell>/...
# which mc_phenotype_assoc_common.sh detects automatically.
set -uo pipefail
. "$(dirname "$0")/mc_phenotype_assoc_common.sh"

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

say "done; results under $out_base"
