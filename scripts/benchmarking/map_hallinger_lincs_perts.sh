#!/bin/bash
# Map the 1,042 Hallinger tier-1 chemicals (structures resolved from names via the
# DSSTox CCD dump by parse_hallinger_supp.py) to the LINCS HepG2 perturbations.

repo="$(cd "$(dirname "$0")/../.." && pwd)"

pert="$repo/data/LINCS/lincs_unique_perts_HEPG2.tsv"
val="$repo/data/validation/hallinger/hallinger_tier1_all1042_annotated.csv"

out_dir="$repo/data/validation/hallinger/mapping"
mkdir -p "$out_dir"

Rscript "$(dirname "$0")/map_smiles.R" \
    --pert "$pert" \
    --val "$val" \
    --pert_smiles_col "smiles" \
    --pert_pertid_col "pert_id" \
    --pert_pertname_col "pertname" \
    --val_smiles_col "SMILES" \
    --val_pertid_col "DTXSID" \
    --val_name_col "chnm" \
    --out_dir "$out_dir"
