#!/bin/bash

# Expression signatures on the MetChange network (networks/combined) with the current recipe, with its
# shuffled-signature control and the comparison. Array job: one fold per task. Seed 0: compare with
# expr_noGlobals_lr1e-3_min15 (same signatures, networks/gene) and metchange_rank_noGlobals_lr1e-3_min15 (same network).
# Why: separates the network from the signature as the cause of MetChange's deficit against expr
# (reports/2026-11-14_seed_tf_expr.md, 10.4).

bash submit_with_shuffle_control.sh expr_combined_noGlobals_lr1e-3_min15
