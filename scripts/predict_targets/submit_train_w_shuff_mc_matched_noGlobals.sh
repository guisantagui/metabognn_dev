#!/bin/bash

# Matched negatives (K = 4, weight 0.5, epochs selected on validation spec) without the per-perturbation inputs
# (--no-globals: relative dose, time, signature size), with its shuffled-signature control and the comparison.
# See reports/2026-11-13_A549_matched_negatives.md, section 5.5 and next steps.

bash submit_with_shuffle_control.sh mc_matched_noGlobals
