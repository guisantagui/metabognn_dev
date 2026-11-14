#!/bin/bash

# The reference recipe with a warm-up before early stopping (--min-epochs 15 --patience 15 --ffn-epochs 80), on
# MetChange (priority), expression and TF (positive control: the published Menacho et al. setting) signatures, each
# with its shuffled-signature control and the comparison. Array jobs: one fold per task.
# Why: about half of the tf / expr FNNs kept epoch 1-2 under patience 10 (validation AUPRC peaks early, dips and
# recovers after the first lr halvings); expr folds 0-1 collapsed (report 03, reports/2026-11-14_seed_tf_expr.md).

bash submit_with_shuffle_control.sh mc_noGlobals_lr1e-3_min15 expr_noGlobals_lr1e-3_min15 tf_noGlobals_lr1e-3_min15
