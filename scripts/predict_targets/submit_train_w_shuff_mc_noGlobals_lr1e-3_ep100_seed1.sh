#!/bin/bash

# Robustness checks of the reference recipe (metchange_rank_noGlobals_lr1e-3), each with its shuffled-signature
# control and the comparison. Array jobs: one fold per task.
# Why (reports/2026-11-14_lr_noglobals_matched.md, 6-7):
#   - ep100: --ffn-epochs 100 instead of 60; one ensemble member kept epoch 53 of 60, so the cap may still bind.
#   - seed1: --seed 1 (new folds, initialisation and permutation); sizes run-to-run variability before building on
#     the result.

bash submit_with_shuffle_control.sh mc_noGlobals_lr1e-3_ep100 mc_noGlobals_lr1e-3_seed1
