#!/bin/bash

proj_dir="/scratch/users/gsantamaria/test_large_files/metabognn_dev/"
cell_file="/home/users/gsantamaria/projects/metabognn_dev/data/LINCS/core_cell_lines.tsv"

mkdir -p logs

while IFS=$'\t' read -r cell; do
    echo "Launching ${cell} sign fix..."
    sbatch \
        --job-name="corr_sign_${cell}" \
        --output="logs/corr_sign_${cell}.out" \
        --error="logs/corr_sign_${cell}.err" \
        launch_fix_limma_reference_sign_expr_cell.sh \
        "$proj_dir" \
        "$cell"
done < "$cell_file"