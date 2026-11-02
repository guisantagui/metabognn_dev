################################################################################
####################----------------------------------------####################
##########----------            MetaboReprogram             ----------##########
# Prepare metchange input matrices: takes one large matrix pulled out from     #
# LINCS GCTX file, converts HGNC symbol IDs to ENSEMBL gene IDs using HGNC's   #
# dictionary, brings expression to linear space (LINCS cycleGAN data is        #
# inferred quantile-normalized log2 counts), removes genes without ENSEMBL     #
# gene IDs and sums eventual duplicates.                                       #
#                                                                              #
################################################################################

import argparse
import math
import sys
import numpy as np
import pandas as pd
from os import makedirs, getcwd
from os.path import join, basename, dirname, abspath
from cobra.io import read_sbml_model

sys.path.insert(0, join(dirname(abspath(__file__)), "..", "parse_lincs_data"))
from gene_id_mapping import symbol_to_ensembl_map, read_matrix, gem_symbol_overrides  # noqa: E402

# Terminal argument parser
################################################################################
parser = argparse.ArgumentParser(
    description='Prepare LINCS CSV files for running MetChange.'
)

parser.add_argument(
        '--expr',
        type=str,
        help='Expression file.'
)

parser.add_argument(
        '--input_type',
        type=str,
        help='Can be expression expression or penalty. If penalty chosen, all the conversion to ensembl, reversal of the log2-transform and rank-based penalty calculation part is skipped.',
        default="expression"
)

parser.add_argument(
        '--mdat',
        type=str,
        help='Metadata file.'
)

parser.add_argument(
        '--gem',
        type=str,
        help='GEM file.',
        default=None
)

parser.add_argument(
        '--hgnc',
        type=str,
        help='HGNC file.',
        default=None
)

parser.add_argument(
        '--batch_size',
        type=int,
        help='Size of the batch, for splitting the CSV into chunks.'
)

parser.add_argument(
        '--n_batches_to_retain',
        type=int,
        help='The number of batches to actually write. Not intended to be used out than for tests/debugging.',
        default=None
)

parser.add_argument('--out_dir', 
	                  type=str,
                    help='Directory where files will be stored.', 
                    default=getcwd())

args = parser.parse_args()

# Directory stuff
################################################################################

expr_f = args.expr
input_type = args.input_type
mdat_f = args.mdat
hgnc_f = args.hgnc
gem_f = args.gem
batch_size = args.batch_size
n_batches_to_retain = args.n_batches_to_retain
out_dir = args.out_dir

makedirs(out_dir, exist_ok=True)

################################################################################
# Load data
################################################################################

# Metadata rows are copied verbatim into the batch files; the id column is
# only used to check that expression columns are in the same order.
print("Loading metadata from %s..." % mdat_f)
with open(mdat_f) as fh:
    mdat_lines = fh.readlines()
mdat_ids = pd.read_csv(mdat_f, dtype=str, keep_default_na=False)["id"].tolist()

if input_type == "expression":
    # Symbols -> ENSEMBL (aliases / previous symbols only if unambiguous), only
    # GEM genes kept while reading, duplicates summed in linear space.
    print("Loading gem %s..." % gem_f)
    gem = read_sbml_model(gem_f)
    gem_genes = set(g.id for g in gem.genes)
    print("Loading HGNC data from %s..." % hgnc_f)
    mapping = symbol_to_ensembl_map(pd.read_csv(hgnc_f, sep='\t', dtype=str), gem_symbol_overrides(gem))
    print("Expression type selected. Translation to ENSEMBL, sum of duplicate IDs, GEM filtering and split in batches will be done.")
    gene_ids, samples, X, report = read_matrix(expr_f, id_col="V1", mapping=mapping, keep=gem_genes, is_log2=True)
    print("%d GEM genes kept (%d rows merged), %d of %d GEM genes missing" % (
        len(gene_ids), (report.status == "merged").sum(), len(gem_genes - set(gene_ids)), len(gem_genes)))
    # Revert log2 expression to count-like (linear) expression
    np.exp2(X, out=X)
elif input_type == "penalty":
    print("Penalty type selected. Only split in batches will be done.")
    gene_ids, samples, X, _ = read_matrix(expr_f, is_log2=False)
else:
    sys.exit("Unsuported input_type.")
print("%d genes x %d samples" % X.shape)

# Expression columns must follow the metadata rows: batches are cut by position.
if samples != mdat_ids:
    missing = sorted(set(mdat_ids) - set(samples))
    extra = sorted(set(samples) - set(mdat_ids))
    if missing or extra:
        sys.exit("Expression columns and metadata ids differ: %d ids only in metadata (e.g. %s), "
                 "%d only in expression (e.g. %s)" % (len(missing), missing[:3], len(extra), extra[:3]))
    print("Expression columns are in a different order than the metadata: reordering them.")
    order = pd.Series(range(len(samples)), index=samples)[mdat_ids].to_numpy()
    X = X[:, order]
    samples = mdat_ids

# Save data
################################################################################

n_samples = len(mdat_ids)
n_batches = math.ceil(n_samples / batch_size)   # the last batch may be smaller
if n_batches_to_retain is not None:
    n_batches = min(n_batches, n_batches_to_retain)

for i in range(n_batches):
    start = i * batch_size
    end = min((i + 1) * batch_size, n_samples)

    mdat_name = join(
        out_dir,
        "%s_batch_%s.csv" % (basename(mdat_f).split('.')[0], i)
    )

    expr_name = join(
        out_dir,
        "%s_batch_%s.csv" % (basename(expr_f).split('.')[0], i)
    )

    print("batch %s: %s to %s" % (i, start, end))

    with open(mdat_name, "w") as fh:
        fh.writelines([mdat_lines[0]] + mdat_lines[1 + start:1 + end])
    pd.DataFrame(X[:, start:end], index=gene_ids, columns=samples[start:end]).to_csv(
        expr_name, index_label="ensembl_gene_id", float_format="%.7g")
