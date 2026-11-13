#!/bin/bash

# Lower learning rate (0.001), more patience (10) and more epochs (60), without the dose/time inputs (--no-globals),
# unmatched and matched, each with its shuffled-signature control and the comparison. Array jobs: one fold per task.
# Why: training stopped after 1-3 epochs in several folds of metchange_rank_noGlobals
# (reports/2026-11-13_A549_matched_negatives.md, 5.5).

bash submit_with_shuffle_control.sh mc_noGlobals_lr1e-3 mc_matched_noGlobals_lr1e-3
