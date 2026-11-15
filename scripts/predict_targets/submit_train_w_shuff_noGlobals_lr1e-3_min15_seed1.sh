#!/bin/bash

# Second seed (--seed 1) of the current recipe (--no-globals --lr 0.001 --min-epochs 15 --patience 15
# --ffn-epochs 80) on MetChange, expression and TF (positive control) signatures, each with its shuffled-signature
# control and the comparison. Array jobs: one fold per task.
# Why: every current-recipe result is one seed; the MetChange-vs-expr difference (3.6 points of per-drug gap) must be
# larger than the seed spread (~1.3 points) (reports/2026-11-14_seed_tf_expr.md, 9.3).

bash submit_with_shuffle_control.sh mc_noGlobals_lr1e-3_min15_seed1 expr_noGlobals_lr1e-3_min15_seed1 \
    tf_noGlobals_lr1e-3_min15_seed1
