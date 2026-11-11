#!/bin/bash

proj_dir="/scratch/users/gsantamaria/test_large_files/metabognn_dev/"

labels_dir="${proj_dir}/results/labels/"

cell_list=("A375" "A549" "HA1E" "HCC515" "HEPG2" "HT29" "MCF7" "VCAP")

for cell in "${cell_list[@]}"; do
    cd ${labels_dir}/${cell}
    mkdir expr/uncorrected
    mv expr/log2FC.csv expr/basal.csv expr/uncorrected/          # keep the originals
    #mv expr/checkpoints expr/uncorrected/ 2>/dev/null            # old per-perturbation limma results (wrong signs)
    rm -r expr/checkpoints
    cp expr_signfixed/log2FC.csv expr_signfixed/basal.csv expr_signfixed/flipped.tsv expr_signfixed/flipped_R.tsv expr/
    echo "reference level: ctrl (sign-fixed by fix_limma_reference_sign.py)" > expr/reference_level.txt
done
