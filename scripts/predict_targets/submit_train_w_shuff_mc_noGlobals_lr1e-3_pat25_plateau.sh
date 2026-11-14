#!/bin/bash

# MetChange, two ways of letting training get past the StepLR halvings, each with its shuffled-signature control and
# the comparison. Array jobs: one fold per task. Same seed (0) and splits as metchange_rank_noGlobals_lr1e-3.
#   - min15_pat25: --patience 25 (StepLR unchanged): every member sees at least two more halvings after its best epoch.
#   - plateau: --lr-schedule plateau --plateau-patience 3 --patience 20: lr reductions tied to validation.
# Why: 29/60 MetChange FNNs kept the first epoch after an lr halving and 52/60 stopped by patience
# (reports/2026-11-14_seed_tf_expr.md, section 8).

bash submit_with_shuffle_control.sh mc_noGlobals_lr1e-3_min15_pat25 mc_noGlobals_lr1e-3_plateau
