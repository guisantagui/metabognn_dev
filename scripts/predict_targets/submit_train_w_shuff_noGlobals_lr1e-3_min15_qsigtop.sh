#!/bin/bash

# The current recipe (_min15, seed 0) with signatures = features with q < 0.05, then the top 100 of them by |log2FC|
# (--qvalue 0.05 --qvalue-top-n --top-n 100). _min15 reads the top 100 by |log2FC| regardless of significance
# (non-significant features enter whenever they have large fold changes); _qsig a random 100 of the significant
# ones, which lost most of the signal (expr: 1.6 vs 7.8 points; reports/2026-11-17_signature_selection_scale.md,
# 5.2-5.3). MetChange (corrected scratch labels), expression and TF, each with its shuffled control and comparison.
# Perturbations with < 3 significant features are dropped (MetChange keeps ~1,791 of 3,970), so compare with
# _min15 on the drugs / perturbations they share.
# MetChange also with --signature-scale sd (_qsigtop_sdscale): the top 100 significant by SD-scaled |log2FC|, so it
# is not a fixed set of high-variance metabolites (report 05, sections 2-4). Compare with _qsigtop (scale effect,
# same perturbations) and _min15_sdscale (significance effect).

bash submit_with_shuffle_control.sh mc_noGlobals_lr1e-3_min15_qsigtop expr_noGlobals_lr1e-3_min15_qsigtop \
    tf_noGlobals_lr1e-3_min15_qsigtop mc_noGlobals_lr1e-3_min15_qsigtop_sdscale
