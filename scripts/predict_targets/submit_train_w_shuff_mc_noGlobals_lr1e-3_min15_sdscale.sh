#!/bin/bash

# MetChange, current recipe (_min15, seed 0, corrected scratch labels), with each metabolite's log2FC divided by its
# SD across perturbations before the top-100 selection (--signature-scale sd; not centred, so signs stay relative to
# the control), plus its shuffled control and comparison. Array jobs: one fold per task.
# Why: unscaled, selection frequency follows metabolite SD (Spearman 0.89), 33 metabolites are in > 50% of
# signatures and unrelated drugs share 22% of their selected features (2.9% after scaling)
# (reports/2026-11-17_signature_selection_scale.md). Compare with the corrected metchange_rank_noGlobals_lr1e-3_min15
# (same seed and folds) and with _qsig (scale-free through significance instead).

bash submit_with_shuffle_control.sh mc_noGlobals_lr1e-3_min15_sdscale
