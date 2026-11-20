#!/bin/bash
# Map the Attene-Ramos compounds (PubChem AID 720637, mitochondrial membrane
# potential in HepG2) onto a LINCS cell line's perturbations by structure.
#
#   cell=HEPG2 bash map_attene_lincs_perts.sh
#
# The cell line is the LINCS side of the mapping, not the assay: the phenotype
# was measured in HepG2 whatever cell line is mapped here. Mapping it onto A549
# therefore gives the CROSS-CELL-LINE CONTROL, not a replication - if the
# association is just as strong with another cell line's labels, it is a
# compound property rather than that cell line's metabolism.
#
# Needs data/LINCS/lincs_unique_perts_<cell>.tsv, written by
# scripts/parse_lincs_data/get_lincs_pert_info.R.
set -uo pipefail

cell="${cell:-HEPG2}"
repo="$(cd "$(dirname "$0")/../.." && pwd)"

pert="${pert:-$repo/data/LINCS/lincs_unique_perts_${cell}.tsv}"
val="${val:-$repo/data/validation/attene_ramos/AID_720637_datatable.csv}"
out_dir="${out_dir:-$repo/data/validation/attene_ramos/mapping_${cell}}"

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
    --val_smiles_col "PUBCHEM_EXT_DATASOURCE_SMILES" \
    --val_pertid_col "PUBCHEM_CID" \
    --out_dir "$out_dir"
