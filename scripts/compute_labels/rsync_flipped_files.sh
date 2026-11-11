#!/bin/bash

eval "$(ssh-agent -s)"
ssh-add --apple-use-keychain ~/.ssh/id_ed25519

labels_dir="/scratch/users/gsantamaria/test_large_files/metabognn_dev/results/labels/"
result_dir="/Users/guillem.santamaria/Documents/postdoc/comput/metabognn_dev/results/labels/"

cell_list="/Users/guillem.santamaria/Documents/postdoc/comput/metabognn_dev/data/LINCS/core_cell_lines.tsv"

while IFS=$'\t' read -r cell; do
    labels_dir_cell="${labels_dir}/${cell}/expr_signfixed/"
    labels_dir_out="${result_dir}/${cell}/expr_signfixed/"
    mkdir -p "$labels_dir_out"
    rsync -avzu iris-cluster:${labels_dir_cell}/flipped* "$labels_dir_out"
done < "$cell_list"
