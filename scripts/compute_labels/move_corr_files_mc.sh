#!/bin/bash

proj_dir="/scratch/users/gsantamaria/test_large_files/metabognn_dev/"

labels_dir="${proj_dir}/results/labels/"

cd ${labels_dir}/A549/metchange
if [[ ! -f rank/reference_level.txt && -s rank_signfixed/log2FC.csv && -s rank_signfixed/basal.csv ]]; then
    mkdir -p rank/uncorrected
    mv rank/log2FC.csv rank/basal.csv rank/uncorrected/
    rm -rf rank/checkpoints
    cp rank_signfixed/{log2FC.csv,basal.csv,flipped.tsv,flipped_R.tsv} rank/
    echo "reference level: ctrl (sign-fixed by fix_limma_reference_sign.py)" > rank/reference_level.txt
fi