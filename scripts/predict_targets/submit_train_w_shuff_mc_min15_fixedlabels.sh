#!/bin/bash

# Retrain MetChange with the current recipe (_min15, seeds 0 and 1), each with its shuffled control and comparison,
# on the corrected labels (scratch metchange/rank: rerun of the fixed limma_DE.R, reference_level.txt "ctrl (explicit)").
# Why: every earlier MetChange run read the home labels, which are the uncorrected limma output (55% of perturbations
# sign-flipped; reports/2026-11-16_labels_network_coverage.md, report 04). The launchers already read scratch labels.
#
# The new runs use the same out_dir names as the buggy ones, so those are first renamed (kept for reference):
# <run> -> <run>_signbug and <run>_shuff -> <run>_signbug_shuff, which keeps the "<real>_shuff" pairing that
# compare_runs.py / check_metabolic_coverage.py expect. Nothing is overwritten: the script stops before moving
# anything if a target already exists.

runs_dir="/scratch/users/gsantamaria/test_large_files/metabognn_dev/results/predict_targets/A549"

moves=()
for run in metchange_rank_noGlobals_lr1e-3_min15 metchange_rank_noGlobals_lr1e-3_min15_seed1; do
    for sfx in "" "_shuff"; do
        src="${runs_dir}/${run}${sfx}"; dst="${runs_dir}/${run}_signbug${sfx}"
        [[ -e "$src" ]] || continue
        [[ -e "$dst" ]] && { echo "${dst} already exists; nothing moved or submitted. Resolve by hand."; exit 1; }
        moves+=("${src}|${dst}")
    done
done
for m in "${moves[@]}"; do
    mv "${m%%|*}" "${m#*|}" && echo "moved ${m%%|*} -> ${m#*|}" || { echo "mv failed"; exit 1; }
done

bash submit_with_shuffle_control.sh mc_noGlobals_lr1e-3_min15 mc_noGlobals_lr1e-3_min15_seed1
