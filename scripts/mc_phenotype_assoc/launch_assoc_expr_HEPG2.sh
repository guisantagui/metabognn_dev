#!/bin/bash
# THE BASELINE: the identical test on the expression labels with the subsystems'
# GPR gene sets. If the transcriptome separates the phenotype groups as well as
# MetChange does, the stoichiometric reading has added nothing - the same
# question the kNN baseline asked of the model (report 06, 6.1).
#
# The GEM states its gene_reaction_rule in Ensembl IDs and the expression label
# matrices are indexed by Ensembl IDs, so no identifier mapping is involved.
#
# Needs gene_sets.tsv:
#   python3 -I scripts/mc_phenotype_assoc/build_metabolite_sets.py --emit genes \
#       --features <expr labels>/log2FC.csv
#
# NOTE the expression label directory must hold log2FC.csv for the value/sign
# modalities; with only qvals.csv present this runs the sig modality alone. It
# also works without fit_info_nsig.csv: drug/dose/time are then parsed from the
# column names and n_sig is counted from the q-values.
set -uo pipefail
. "$(dirname "$0")/mc_phenotype_assoc_common.sh"

# see the note in launch_assoc_metchange_HEPG2.sh: sig runs on the responsive
# subset only (--min-nsig 1), value/sign on everything.
sig_modes="${sig_modes:---mode sig --min-nsig 1}"
if [ -f "$expr_labels/log2FC.csv" ]; then
    modes="${modes:---mode value --mode sign}"
else
    echo "no log2FC.csv in $expr_labels - the value/sign modalities are skipped" >&2
    modes=""
fi

if [ -n "$modes" ]; then
    run_assoc attene_expr "$expr_labels" "$gene_sets" \
        "$pheno_dir/attene_ramos_mmp.tsv" MMP_down inactive \
        $modes --extra-contrast MMP_up:inactive --extra-contrast cytotox:inactive
    run_strat attene_expr MMP_down
    run_strat attene_expr MMP_up

    run_assoc hallinger_expr "$expr_labels" "$gene_sets" \
        "$pheno_dir/hallinger_ocr.tsv" active inactive $modes
    run_strat hallinger_expr active
fi

run_assoc attene_expr_sig "$expr_labels" "$gene_sets" \
    "$pheno_dir/attene_ramos_mmp.tsv" MMP_down inactive \
    $sig_modes --extra-contrast MMP_up:inactive
run_assoc hallinger_expr_sig "$expr_labels" "$gene_sets" \
    "$pheno_dir/hallinger_ocr.tsv" active inactive $sig_modes

say "done; compare set_results.tsv here with the metchange runs set by set"
