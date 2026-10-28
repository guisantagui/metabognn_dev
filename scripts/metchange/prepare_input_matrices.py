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
import pandas as pd
import pyarrow
from os import makedirs, getcwd
from os.path import join, basename
import polars as pl
from collections import defaultdict
import cobra
from cobra.io import read_sbml_model

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

# Load data
################################################################################

print("Loading metadata from %s..." % mdat_f)
mdat = pl.read_csv(mdat_f)
print("Loading expression data from %s..." % expr_f)
expr = pl.read_csv(expr_f)

n_batches = round(mdat.shape[0] / batch_size)

# Parse data
################################################################################

if input_type == "expression":
    # If expression, these files will be needed too.
    print("Loading gem %s..." % gem_f)
    gem = read_sbml_model(gem_f)
    print("Loading HGNC data from %s..." % hgnc_f)
    hgnc = pd.read_csv(hgnc_f, sep = '\t')

    print("Expression type selected. Translation to ENSEMBL, sum of duplicate IDs, GEM filtering and split in batches will be done.")
    # Create symbol to ENSEMBL dictionary
    symbol_to_ensembl = (
        hgnc.groupby("symbol")["ensembl_gene_id"]
        .apply(lambda x: "|".join(x.dropna().astype(str)))
        .to_dict()
    )

    alias_to_ensembl = defaultdict(list)

    for _, row in hgnc.dropna(subset=["alias_symbol"]).iterrows():
        for alias in row["alias_symbol"].split("|"):
            alias_to_ensembl[alias].append(str(row["ensembl_gene_id"]))

    alias_to_ensembl = {
        k: "|".join(set(v)) for k, v in alias_to_ensembl.items()
    }

    ensembl_list = []

    for gene in expr['V1'].to_list():
        if gene in symbol_to_ensembl:
            ensembl_list.append(symbol_to_ensembl[gene])
        else:
            ensembl_list.append(alias_to_ensembl.get(gene, ""))
    symb_to_ensembl = pd.DataFrame(
        {
            'hgnc_symbol': expr['V1'].to_list(),
            'ensembl_gene_id': ensembl_list
        }
    )

    # Revert log2 expression to count-like expression, add ensembl column,
    # remove genes without mapped ensembl and sum duplicated ENSEMBL IDs
    expr = expr.with_columns(
        pl.exclude("V1").pow(2)
    )

    expr = expr.with_columns(
        pl.Series("ensembl_gene_id", ensembl_list)
    )

    expr = expr.filter(
        pl.col("ensembl_gene_id") != ""
    )

    expr = expr.group_by("ensembl_gene_id").agg(
        pl.exclude(["V1", "ensembl_gene_id"]).sum()
    )

    gem_genes = [x.id for x in gem.genes]

    expr = expr.filter(
        pl.col("ensembl_gene_id").is_in(gem_genes)
    )
elif input_type == "penalty":
    print("Penalty type selected. Only split in batches will be done.")
else:
    sys.exit("Unsuported input_type.")

# Save data
################################################################################

gene_ids = expr.select("ensembl_gene_id")

if n_batches_to_retain != None:
    end_loop = n_batches_to_retain + 1
else:
    end_loop = n_batches

for i in range(0, end_loop):
    start = i * batch_size
    end = min((i + 1) * batch_size, mdat.shape[0])

    mdat_batch = mdat[start:end, :]

    # Skip ensembl_gene_id when slicing expression samples
    expr_batch = expr[:, (start + 1):(end + 1)]

    # Add gene IDs back as the first column
    expr_batch = gene_ids.hstack(expr_batch)

    mdat_name = join(
        out_dir,
        "%s_batch_%s.csv" % (basename(mdat_f).split('.')[0], i)
    )

    expr_name = join(
        out_dir,
        "%s_batch_%s.csv" % (basename(expr_f).split('.')[0], i)
    )

    print("batch %s: %s to %s" % (i, start, end))

    mdat_batch.write_csv(mdat_name)
    expr_batch.write_csv(expr_name)