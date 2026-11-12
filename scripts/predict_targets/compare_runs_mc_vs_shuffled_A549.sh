#!/bin/bash

cell_dir="/home/users/gsantamaria/projects/metabognn_dev/results/predict_targets/A549/"
real="${cell_dir}/metchange_rank/"
lin_shuf="${cell_dir}/metchange_rank_linOnly_shuff"
shuf="${cell_dir}/metchange_rank_shuff/"
out_dir="${cell_dir}/mc_real_vs_suff/"

mkdir -p "$out_dir"

sbatch launch_compare_runs.sh \
    --run real=$real \
    --run lin_shuff="$lin_shuf" \
    --run shuffled=$shuf \
    --control shuffled \
    --out-dir "$out_dir"