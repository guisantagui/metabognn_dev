#!/bin/bash
# Build the inputs the association launchers need, for one cell line.
#
#   cell=HEPG2 bash prepare_assoc_inputs.sh    # phenotypes for one cell line
#   cell=ALL   bash prepare_assoc_inputs.sh    # ONE union table for every cell line
#
# cell=ALL maps the external datasets against lincs_unique_perts.tsv (the union
# over cell lines) and writes <assoc_base>/phenotypes_ALL/. The association
# script keeps only the drugs that have labels in whichever cell line it runs
# on, so that one table covers all of them and needs a single map_smiles.R run
# instead of one per cell line. A per-cell directory, when present, takes
# precedence over it.
#
# Writes:
#   <assoc_base>/metabolite_sets/{metabolite_sets,gene_sets,metabolite_sets_negctrl}.tsv
#       PROJECT-level: built from Human-GEM alone, so cell-independent. They are
#       deliberately built WITHOUT --features, i.e. not restricted to one label
#       matrix; mc_phenotype_assoc.py drops members absent from whatever matrix
#       it is given and reports how many it found. (Restricting them to a cell
#       line is what used to make them cell-coupled.)
#   <assoc_base>/<cell>/phenotypes/{attene_ramos_mmp,hallinger_ocr,hallinger_ocr_bymech}.tsv
#       PER-CELL: keyed on this cell line's LINCS perturbation names.
#
# The phenotype tables need a structure mapping of the external datasets onto
# this cell line's perturbations, which needs R + rcdk:
#   cell=<cell> bash scripts/benchmarking/map_attene_lincs_perts.sh
#   cell=<cell> bash scripts/benchmarking/map_hallinger_lincs_perts.sh
# and those in turn need data/LINCS/lincs_unique_perts_<cell>.tsv from
# scripts/parse_lincs_data/get_lincs_pert_info.R. This script says so and stops
# if the mapping is absent rather than falling back to another cell line's.
set -uo pipefail
. "$(dirname "$0")/mc_phenotype_assoc_common.sh"
provenance
require_python || exit 1

mkdir -p "$sets_dir" "$pheno_dir"

# ---------------------------------------------------------------- sets (shared)
if [ "${rebuild_sets:-false}" = "true" ] || [ ! -f "$met_sets" ]; then
    say "metabolite sets -> $sets_dir"
    "$py" -I "$here_dir/build_metabolite_sets.py" --out-dir "$sets_dir"
    say "declared negative-control sets -> $sets_dir"
    "$py" -I "$here_dir/build_metabolite_sets.py" --emit negative-controls --out-dir "$sets_dir"
    say "GPR gene sets (expression baseline) -> $sets_dir"
    "$py" -I "$here_dir/build_metabolite_sets.py" --emit genes --out-dir "$sets_dir"
else
    echo "metabolite sets already present ($met_sets); rebuild_sets=true to force" >&2
fi

# ------------------------------------------------------------ phenotypes (cell)
say "phenotype tables for cell=$cell -> $pheno_dir"
if ! "$py" -I "$here_dir/build_phenotype_tables.py" --cell "$cell" --out-dir "$pheno_dir"; then
    cat <<MSG >&2

Phenotype tables could not be built for cell=$cell. Most likely the structure
mapping for this cell line does not exist yet. Run, in order:

    Rscript scripts/parse_lincs_data/get_lincs_pert_info.R   # -> lincs_unique_perts_${cell}.tsv
    cell=$cell bash scripts/benchmarking/map_attene_lincs_perts.sh
    cell=$cell bash scripts/benchmarking/map_hallinger_lincs_perts.sh
    cell=$cell bash $(basename "$0")

NOTE the phenotypes are measured in HepG2 whatever cell line is mapped. Mapping
them onto another cell line is the CROSS-CELL-LINE CONTROL, not a replication.
MSG
    exit 1
fi

# by-mechanism split of the Hallinger table (ETCi / uncoupler / ATP synthase)
"$py" -I - "$pheno_dir" <<'PYEOF'
import csv, os, sys, collections
d = sys.argv[1]
src = os.path.join(d, 'hallinger_ocr.tsv')
if os.path.exists(src):
    rows = list(csv.DictReader(open(src), delimiter='\t'))
    for r in rows:
        if r['group'] == 'active':
            m = r['subgroup']
            r['group'] = {'ETCi': 'ETCi', 'Uncoupler': 'uncoupler',
                          'ATP Synthase': 'atpsynth'}.get(m, 'other_active')
    out = os.path.join(d, 'hallinger_ocr_bymech.tsv')
    with open(out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=rows[0].keys(), delimiter='\t')
        w.writeheader(); w.writerows(rows)
    print('  hallinger_ocr_bymech.tsv: %s' % dict(collections.Counter(r['group'] for r in rows)))
PYEOF

say "done. sets: $sets_dir ; phenotypes: $pheno_dir"
