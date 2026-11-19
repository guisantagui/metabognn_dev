#!/bin/bash

pert="/Users/guillem.santamaria/Documents/postdoc/comput/metabognn_dev/data/LINCS/lincs_unique_perts.tsv"

val="/Users/guillem.santamaria/Documents/postdoc/comput/metabognn_dev/data/validation/attene_ramos/AID_720637_datatable.csv"

val_smiles_col="PUBCHEM_EXT_DATASOURCE_SMILES"
pert_smiles_col="smiles"
pert_pertid_col="pert_id"
pert_pertname_col="pertname"
val_pertid_col="PUBCHEM_CID"

out_dir="$(dirname "$val")"

Rscript map_smiles.R \
    --pert "$pert" \
    --val "$val" \
    --pert_smiles_col "$pert_smiles_col" \
    --pert_pertid_col "$pert_pertid_col" \
    --pert_pertname_col "$pert_pertname_col" \
    --val_smiles_col "$val_smiles_col" \
    --val_pertid_col "$val_pertid_col" \
    --out_dir "$out_dir"

