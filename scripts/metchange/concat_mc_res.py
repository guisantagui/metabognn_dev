################################################################################
# MetaboPert: Concatenate MetChange batch matrices, compute deltas by          #
# subtracting the mean of the matched controls inside of each batch, obtain    #
# the average per perturbation and save result and metadata in the output      #
# directory.                                                                   #
################################################################################

import os
import polars as pl
from pathlib import Path
import re
import numpy as np
import argparse

# Terminal argument parser
################################################################################
parser = argparse.ArgumentParser(
    description='Concatenate MetChange batch matrices, compute deltas by subtracting the mean of the matched controls inside of each batch, obtain the average per perturbation and save result and metadata in the output directory.'
)

parser.add_argument(
        '--mc_dir',
        type=str,
        help='Directory with metchange batch files.'
)

parser.add_argument(
        '--out_dir', 
	type=str,
        help='Directory where files will be stored.', 
        default=os.getcwd()
)

args = parser.parse_args()

# Directory stuff
################################################################################

mc_dir = args.mc_dir
out_dir = args.out_dir

os.makedirs(out_dir, exist_ok=True)

# Data loading
################################################################################

# Get metchange files, filtering out any concatenated file or the .npy of LP1
mc_files = [
        x for x in os.listdir(mc_dir) if ".npy" not in x and ("_concat.csv" not in x) and ("_delta.csv" not in x)
]

mc_files = [x for x in mc_files if "rank_penalties.csv" not in x]

# Concatenate the preocessing batches metChange matrices into a single one
mc_all = None
for i in range(0, len(mc_files)):
    mc = re.sub(r'(?<=batch_)\d+', str(i), mc_files[0])
    print("processing %s..." % mc)
    mc_df = pl.read_csv(os.path.join(mc_dir, mc))
    batch = mc.split('_')[4]
    pattern = rf"_{batch}\."
    mc_df = mc_df.rename({mc_df.columns[0]: "compound"})
    if mc_all is None:
        mc_all = mc_df
    else:
        mc_all = mc_all.hstack(mc_df.drop("compound"))

# Remove compounds where all the samples are None for a metabolite. These
# are the metabolites that can't be produced by the model at all under the
# specified media conditions.
all_none = mc_all.select(
    pl.all_horizontal(pl.exclude("compound").is_null()).alias("all_none")
)["all_none"].to_list()
mc_all = mc_all.filter(~pl.Series(all_none))

# Save dataframe
################################################################################

outname_concat = re.sub(
    r'_batch_\d+',
    '',
    mc_files[0]
).replace('.csv', '_concat.csv')

outname_delta = outname_concat.replace('.csv', '_delta.csv')

mc_all.write_csv(os.path.join(out_dir, outname_concat))
print("%s saved at %s." % (outname_concat, out_dir))