################################################################################
####################----------------------------------------####################
##########----------            MetaboReprogram             ----------##########
# Convert gene symbols to ensembl and aggregate duplicated. LINCS database     #
# comes with SYMBOL IDs, and Human2-GEM has ENSEMBL in GPRs. So, conversion is #
# required.                                                                    #
#                                                                              #
# The matrix (genes x samples) is streamed in row chunks, so memory stays at a #
# few GB whatever the number of samples: only the rows of genes that several   #
# symbols map to are held until their group is complete.                       #
################################################################################

import os
import argparse
import csv
from collections import defaultdict

import numpy as np
import pandas as pd

from gene_id_mapping import symbol_to_ensembl_map, merge_rows, gem_symbol_overrides


################################################################################
# Terminal argument parser                                                     #
################################################################################
parser = argparse.ArgumentParser(
    description='Convert the gene symbols of an expression matrix (genes x samples) to ENSEMBL.'
)

parser.add_argument(
        '--expr',
        type=str,
        help='Expression file.'
)

parser.add_argument(
        '--hgnc',
        type=str,
        help='HGNC file, for translating symbols to ENSEMBL.'
)

parser.add_argument(
        '--symbol_col',
        type=str,
        help='Column with the gene symbol IDs.'
)

parser.add_argument(
        '--is_log2',
        action="store_true",
        help='If input expression data is log transformed.'
)

parser.add_argument(
        '--gem',
        type=str,
        default=None,
        help='GEM file (optional). Its genes keep their own ENSEMBL ids where HGNC lists another one.'
)

parser.add_argument(
        '--agg',
        choices=["sum", "mean"],
        default="sum",
        help='How rows mapping to the same ENSEMBL id are combined (in linear space if --is_log2).'
)

parser.add_argument(
        '--chunk_rows',
        type=int,
        default=100,
        help='Genes read per chunk (memory ~ chunk_rows x n_samples x 8 bytes).'
)

parser.add_argument(
        '--delete_original',
        action="store_true",
        help='Whether original file should be deleted after ID conversion.'
)

parser.add_argument('--out_dir',
                    type=str,
                    help='Directory where the result will be stored.',
                    default=None)

args = parser.parse_args()

################################################################################
# DIRECTORY STUFF                                                              #
################################################################################
expr_f = args.expr
symbol_col = args.symbol_col
out_dir = args.out_dir if args.out_dir is not None else os.path.dirname(expr_f)
out_name = "%s_ensembl.csv" % os.path.basename(expr_f).split('.')[0]
out_f = os.path.join(out_dir, out_name)

################################################################################
# MAPPING (first pass: symbol column only)                                     #
################################################################################
hgnc = pd.read_csv(args.hgnc, sep='\t', dtype=str)
overrides = None
if args.gem:
    from cobra.io import read_sbml_model
    overrides = gem_symbol_overrides(read_sbml_model(args.gem))
mapping = symbol_to_ensembl_map(hgnc, overrides)

with open(expr_f) as fh:
    header = next(csv.reader(fh))
sym_idx = header.index(symbol_col)
symbols = pd.read_csv(expr_f, usecols=[sym_idx], dtype=str).iloc[:, 0].tolist()
ensembl = [mapping.get(s, "") for s in symbols]
n_rows_per_id = defaultdict(int)
for e in ensembl:
    if e:
        n_rows_per_id[e] += 1
report = pd.DataFrame({
    "symbol": symbols, "ensembl_gene_id": ensembl,
    "status": ["dropped" if not e else "merged" if n_rows_per_id[e] > 1 else "kept" for e in ensembl],
})
report.to_csv(os.path.join(out_dir, "%s_ensembl_mapping.tsv" % os.path.basename(expr_f).split('.')[0]),
              sep="\t", index=False)
print("%d genes: %d kept, %d merged into %d ENSEMBL ids, %d dropped (no unambiguous ENSEMBL id)" % (
    len(symbols), (report.status == "kept").sum(), (report.status == "merged").sum(),
    sum(1 for v in n_rows_per_id.values() if v > 1), (report.status == "dropped").sum()))

################################################################################
# CONVERSION AND SAVING (second pass: row chunks)                              #
################################################################################
samples = [h for i, h in enumerate(header) if i != sym_idx]
pending = defaultdict(list)  # ENSEMBL id -> rows still waiting for the rest of their group
written = 0
print("Converting %s (%d samples) in chunks of %d genes..." % (expr_f, len(samples), args.chunk_rows))
with open(out_f, "w", newline="") as out:
    csv.writer(out).writerow(["ensembl_gene_id"] + samples)
    reader = pd.read_csv(expr_f, chunksize=args.chunk_rows, index_col=sym_idx, dtype={symbol_col: str})
    row = 0
    for chunk in reader:
        values = chunk[samples].to_numpy(dtype=np.float64)
        ids, rows_out = [], []
        for k in range(len(chunk)):
            e = ensembl[row]
            row += 1
            if not e:
                continue
            if n_rows_per_id[e] == 1:
                vals = values[k]
            else:
                pending[e].append(values[k])
                if len(pending[e]) < n_rows_per_id[e]:
                    continue
                vals = merge_rows(pending.pop(e), args.is_log2, args.agg)
            ids.append(e)
            rows_out.append(vals)
        if ids:  # vectorised writer: formatting value by value in Python would take hours
            pd.DataFrame(np.vstack(rows_out), index=ids).to_csv(out, header=False, float_format="%.7g")
            written += len(ids)
        print("  %d / %d genes read" % (row, len(symbols)), flush=True)
assert not pending, "incomplete duplicate groups: %s" % list(pending)[:5]
print("%d ENSEMBL genes written to %s" % (written, out_f))

if args.delete_original:
    print("Deleting original input file...")
    os.remove(expr_f)
    print("%s deleted." % expr_f)
