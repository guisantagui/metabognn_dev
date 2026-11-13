#!/bin/bash

# Two ablations, each with its shuffled-signature control and the comparison against it:
# - mc_noGlobals: unmatched run without the per-perturbation features (relative dose, time, signature size), to test
#   whether the unmatched model's specificity is dose/time (its shuffled control should then sit at 0.5 unstratified).
# - mc_matched_auprc: matched negatives with --select-metric auprc (instead of spec), to separate the effect of the
#   sampling from that of the epoch selection when comparing with metchange_rank and metchange_rank_matched.

bash submit_with_shuffle_control.sh mc_noGlobals mc_matched_auprc
