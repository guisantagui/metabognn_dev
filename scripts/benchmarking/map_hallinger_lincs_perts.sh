#!/bin/bash
# Map the 1,042 Hallinger tier-1 chemicals (structures resolved from names via the
# DSSTox CCD dump by parse_hallinger_supp.py) to the LINCS HepG2 perturbations.

#   cell=HEPG2 bash map_hallinger_lincs_perts.sh
# Needs data/LINCS/lincs_unique_perts_<cell>.tsv (get_lincs_pert_info.R).
set -uo pipefail

cell="${cell:-HEPG2}"
repo="$(cd "$(dirname "$0")/../.." && pwd)"

pert="${pert:-$repo/data/LINCS/lincs_unique_perts_${cell}.tsv}"
val="${val:-$repo/data/validation/hallinger/hallinger_tier1_all1042_annotated.csv}"

out_dir="${out_dir:-$repo/data/validation/hallinger/mapping_${cell}}"
if [ ! -f "$pert" ]; then
    echo "no perturbation table for cell=$cell: $pert" >&2
    echo "generate it with scripts/parse_lincs_data/get_lincs_pert_info.R" >&2
    exit 1
fi
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
