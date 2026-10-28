################################################################################
####################----------------------------------------####################
##########----------            MetaboReprogram             ----------##########
# Convert gene symbols to ensembl and aggregate duplicated. LINCS database     #
# comes with SYMBOL IDs, and Human2-GEM has ENSEMBL in GPRs. So, conversion is #
# required.                                                                    #
#                                                                              #
################################################################################

import os
import polars as pl
import pandas as pd
from collections import defaultdict
import argparse


################################################################################
# MODULE LOADING                                                               #
################################################################################

def convert_ids_ensembl(df, hgnc, symbol_col="V1", is_log2=True):
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
        k: "|".join(sorted(set(v)))
        for k, v in alias_to_ensembl.items()
    }
    ensembl_list = []
    for gene in df[symbol_col].to_list():
        if gene in symbol_to_ensembl:
            ensembl_list.append(symbol_to_ensembl[gene])
        else:
            ensembl_list.append(alias_to_ensembl.get(gene, ""))
    symb_to_ensembl = pd.DataFrame(
        {
            'hgnc_symbol': df[symbol_col].to_list(),
            'ensembl_gene_id': ensembl_list
        }
    )
    # Add ensembl column
    df = df.with_columns(
        pl.Series("ensembl_gene_id", ensembl_list)
    )
    # Remove unmatched rows
    df = df.filter(
        pl.col("ensembl_gene_id") != ""
    )
    if is_log2:
        df = df.with_columns(
            pl.exclude(["ensembl_gene_id", symbol_col]).pow(2)
        )
    # Sum eventual duplicates
    df = (
        df.group_by("ensembl_gene_id")
        .agg(pl.exclude([symbol_col, "ensembl_gene_id"]).sum())
        .sort("ensembl_gene_id")
    )
    if is_log2:
        df = df.with_columns(
            pl.exclude("ensembl_gene_id").sqrt()
        )
    return (df)

################################################################################
# Terminal argument parser                                                     #
################################################################################
parser = argparse.ArgumentParser(
    description='Obtain penalty matrices.'
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
hgnc_f = args.hgnc
symbol_col = args.symbol_col
is_log2 = args.is_log2
delete_original = args.delete_original
out_dir = args.out_dir

if out_dir == None:
    out_dir = os.dirname(expr_f)


################################################################################
# DATA LOADING.                                                                #
################################################################################
print("Loading expression data from %s..." % expr_f)
expr = pl.read_csv(expr_f)
hgnc = pd.read_csv(hgnc_f, sep = '\t')

################################################################################
# CONVERSION AND SAVING                                                        #
################################################################################

print("Converting gene symbols to ENSEMBL...")
expr = convert_ids_ensembl(expr, hgnc, symbol_col=symbol_col, is_log2=is_log2)

out_name = "%s_ensembl.csv" % os.path.basename(expr_f).split('.')[0]

expr.write_csv(os.path.join(out_dir, out_name))

if delete_original:
    print("Deleting original input file...")
    os.remove(expr_f)
    print("%s deleted." % expr_f)