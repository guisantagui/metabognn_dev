#!/bin/bash

# The plateau schedule (--lr-schedule plateau --plateau-patience 3 --patience 20 --min-epochs 15 --ffn-epochs 80) on
# expression and TF (positive control) signatures, each with its shuffled-signature control and the comparison.
# Array jobs: one fold per task. Compare with the _min15 runs (StepLR, same warm-up) and with MetChange _plateau
# (reports/2026-11-14_seed_tf_expr.md, section 8).

bash submit_with_shuffle_control.sh expr_noGlobals_lr1e-3_plateau tf_noGlobals_lr1e-3_plateau
