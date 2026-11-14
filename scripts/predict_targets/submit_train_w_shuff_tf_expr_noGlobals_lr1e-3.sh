#!/bin/bash

# The reference recipe (--no-globals --lr 0.001 --patience 10 --ffn-epochs 60) on the TF-activity and expression
# signatures, each with its shuffled-signature control and the comparison. Array jobs: one fold per task.
# Why: first interpretable comparison across signature types, each run against its own control
# (reports/2026-11-14_lr_noglobals_matched.md, 7.3). Compare with metchange_rank_noGlobals_lr1e-3.

bash submit_with_shuffle_control.sh tf_noGlobals_lr1e-3 expr_noGlobals_lr1e-3
