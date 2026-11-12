#!/bin/bash
cell_file="/home/users/gsantamaria/projects/metabognn_dev/data/LINCS/core_cell_lines.tsv"
proj_dir="/scratch/users/gsantamaria/test_large_files/metabognn_dev/"
workers=64
mode="mc"
penalty="rank"
medium="/home/users/gsantamaria/projects/metabo_pert/data/GEMs/gem_mediums/Human2-GEM/rpmi_medium.json"

mkdir -p "logs/"

while IFS=$'\t' read -r cell; do
    job_name="limma_${cell}_${mode}"
    sbatch \
        --job-name="$job_name" \
        --output="logs/${job_name}.out" \
        --error="logs/${job_name}.err" \
        launch_limma_DE_generic_mc.sh \
            "$proj_dir" \
            "$cell" \
            "$workers" \
            "$mode" \
            "$penalty" \
            "$medium"
done < "$cell_file"