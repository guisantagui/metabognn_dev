

import sys
sys.path.insert(0, 'scripts/parse_lincs_data')

import xlrd
import pandas as pd
from cobra.io import read_sbml_model
from gene_id_mapping import symbol_to_ensembl_map, gem_symbol_overrides
import os
import numpy as np

def isnan(num):
    return num != num

gem = read_sbml_model('data/GEMs/Human-GEM/model/Human-GEM.xml')
nci60_f = "data/validation/human1_nci60/output 2/RNA__RNA_seq_composite_expression.xls"

out_dir="data/validation/human1_nci60/"
mapping = symbol_to_ensembl_map(
    pd.read_csv('data/utilities/hgnc_complete_set.txt', sep='\t', dtype=str),
    gem_symbol_overrides(gem),
)

nci60 = pd.read_excel(nci60_f, engine="xlrd")

i = 0
while isnan(nci60.iloc[i, 3]):
    i += 1

cols = nci60.iloc[i, :].to_list()

nci60 = nci60.iloc[(i + 1):nci60.shape[0], :]
nci60.columns = cols

nci60_ensembl = nci60['Gene name d'].map(mapping)

print(
    "%f %% mapped to ENSEMBL gene IDs." % (
        sum([not isnan(x) for x in nci60_ensembl]) / nci60.shape[0] * 100
    )
)

nci60 = nci60.iloc[:, 6:nci60.shape[1]]

nci60.insert(loc=0, column="ensembl_gene_id", value=nci60_ensembl)

nci60 = nci60[[not isnan(x) for x in nci60["ensembl_gene_id"]]]

nci60 = nci60.groupby("ensembl_gene_id").max()

nci60.to_csv(os.path.join(out_dir, "nci60_expr_ensembl.csv"), index_label=None)