#!/bin/bash

# The current recipe (_min15: --no-globals --lr 0.001 --min-epochs 15 --patience 15 --ffn-epochs 80, seed 0) with
# significance-thresholded signatures: features with q < 0.05, randomly capped at 100 (--qvalue 0.05
# --max-sources 100), instead of the top 100 by |log2FC| regardless of significance (--ignore-qvalue --top-n 100).
# This is the published Menacho et al. setting (significant DE TFs, max 100). MetChange (corrected scratch labels),
# expression and TF, each with its shuffled control and comparison. Array jobs: one fold per task.
# Why: 57% of active perturbations have no significant MetChange feature among the top 100 read, so their signature
# is noise (reports/2026-11-16_labels_network_coverage.md, section 3). Perturbations with < 3 significant features
# are dropped (--min-sources 3), so each signature type is scored on its own perturbation set: compare across
# signature types on the drugs / perturbations they share, and against the _min15 runs of the same signature.

bash submit_with_shuffle_control.sh mc_noGlobals_lr1e-3_min15_qsig expr_noGlobals_lr1e-3_min15_qsig \
    tf_noGlobals_lr1e-3_min15_qsig
