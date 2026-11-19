#!/bin/bash
# Re-map the Schuhknecht drugs to LINCS with the fixed map_smiles.R, to check the
# original mapping (scripts/schuhknecht_map, rcdk "Canonical" + intersect/match:
# no salt stripping, stereochemistry silently dropped, first hit only, NA able to
# match NA). Sheet 2 of the supplement holds "Drug Name" / "Drug SMILES".

repo="$(cd "$(dirname "$0")/../.." && pwd)"
here="$(cd "$(dirname "$0")" && pwd)"

xlsx="${SCHUH_XLSX:-$repo/data/validation/schuhknecht/NBT-RA62370A-Z_NonPooledFoldChanges.xlsx}"
pert="${LINCS_CPDINFO:-/Users/guillem.santamaria/Documents/postdoc/comput/metabo_pert/data/LINCS/compoundinfo_beta.txt}"
out_dir="$repo/data/validation/schuhknecht/mapping_remap"
mkdir -p "$out_dir"

val_csv="$out_dir/schuhknecht_drug_info.csv"
Rscript -e "
  d <- readxl::read_xlsx('$xlsx', sheet = 2)
  d <- d[, c('Drug Name', 'Drug SMILES')]
  names(d) <- c('drug_name', 'drug_smiles')
  write.csv(d, '$val_csv', row.names = FALSE)
  cat(sprintf('schuhknecht sheet 2: %d rows, %d unique names\n', nrow(d), length(unique(d\$drug_name))))
"

Rscript "$here/map_smiles.R" \
    --pert "$pert" \
    --val "$val_csv" \
    --pert_smiles_col "canonical_smiles" \
    --pert_pertid_col "pert_id" \
    --pert_pertname_col "cmap_name" \
    --val_smiles_col "drug_smiles" \
    --val_pertid_col "drug_name" \
    --val_name_col "drug_name" \
    --out_dir "$out_dir"
