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

# Functions
################################################################################
def filt_cols(mc_df, mdat_df, value, mdat_var):
    cols = [
        col for col, p in zip(mc_df.columns[1:], mdat_df[mdat_var])
        if p == value
    ]
    return cols

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
        '--mdat',
        type=str,
        help='Metadata file.'
)

parser.add_argument(
        '--batch_col',
        type=str,
        help='Metadata column with the batch information.',
        default="batch"
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
#mc_dir = "/scratch/users/gsantamaria/test_large_files/metabo_perts/results/metChange/A549_mixPnlt/"
#mdat_f = "/home/users/gsantamaria/projects/metabo_pert/results/gene_penalties/mdat_penalties.csv"
#batch_col = "batch"
#out_dir = "/home/users/gsantamaria/projects/metabo_pert/results/mc_deltas/A549/mixPnlty/"


mc_dir = args.mc_dir
mdat_f = args.mdat
batch_col = args.batch_col
out_dir = args.out_dir

os.makedirs(out_dir, exist_ok=True)

# Data loading
################################################################################

# Load metadata
mdat = pl.read_csv(mdat_f)

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
    #matches = [f for f in mdat_files if re.search(pattern, f)]
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

# Delta calculation
################################################################################
plates = np.unique(list(mdat[batch_col]))
delta_all = None
mdat_delta_all = None
for i in range(0, len(plates)):
    plate = plates[i]
    cols = ["compound"] + filt_cols(mc_all, mdat, plate, batch_col)
    mc_plate = mc_all.select(cols)
    mdat_plate = mdat.filter([x == plate for x in mdat[batch_col]])
    print("Processing batch %i: %s..." % (i, plate))
    for timepoint in np.unique(list(mdat_plate["timepoint"])):
        cols_timepoint = filt_cols(mc_plate, mdat_plate, timepoint, "timepoint")
        mc_plate_tp = mc_plate.select(["compound"] + cols_timepoint)
        mdat_plate_tp = mdat_plate.filter([x == timepoint for x in mdat_plate["timepoint"]])
        ctrl_cols = filt_cols(mc_plate_tp, mdat_plate_tp, "ctrl", "exper_group")
        if len(ctrl_cols) > 0:
            non_ctrl_cols = [
                col for col in mc_plate_tp.columns[1:]
                if col not in ctrl_cols
            ]
            if len(non_ctrl_cols) > 0:
                delta_df = mc_plate_tp.select(
                "compound",
                    *[
                        (pl.col(col) - pl.concat_list(ctrl_cols).list.mean())
                        for col in non_ctrl_cols
                    ]
                )
                if delta_all is None:
                    delta_all = delta_df
                else:
                    delta_all = delta_all.hstack(delta_df.drop("compound"))
                
                order = pl.DataFrame({
                    "id": non_ctrl_cols,
                    "order": range(len(non_ctrl_cols))
                })
                mdat_delta_plate_tp = (
                    mdat_plate_tp
                    .join(order, on="id", how="inner")
                    .sort("order")
                    .drop("order")
                )
                if mdat_delta_all is None:
                    mdat_delta_all = mdat_delta_plate_tp
                else:
                    mdat_delta_all = mdat_delta_all.vstack(mdat_delta_plate_tp)
        else:
            print("No controls found in batch %s. Skipping..." % plate)

# Get average of deltas
################################################################################

print("Computing average per perturbation...")
perturbations = mdat_delta_all["perturbation"].to_list()
sample_cols = delta_all.columns[1:delta_all.shape[1]]

X = delta_all.select(sample_cols).to_numpy()

unique_perts = list(dict.fromkeys(perturbations))

delta_avg = pl.DataFrame({
    "compound": delta_all["compound"],
    **{
        pert: X[:, np.array(perturbations) == pert].mean(axis=1)
        for pert in unique_perts
    }
})

mc_mdat_avg = (
    mdat_delta_all
    .select(["dose", "pertname", "timepoint", "perturbation"])
    .unique()
    .with_columns(
        pl.col("perturbation")
        .replace_strict(
            {p: i for i, p in enumerate(delta_avg.columns)},
            default=None
        )
        .alias("_order")
    )
    .sort("_order")
    .drop("_order")
)

# Save dataframes
################################################################################

outname_concat = re.sub(
    r'_batch_\d+',
    '',
    mc_files[0]
).replace('.csv', '_concat.csv')

outname_delta = outname_concat.replace('.csv', '_delta.csv')
outname_mdat_delta = os.path.basename(mdat_f).replace(".csv", "_delta.csv")
outname_delta_avg = outname_concat.replace('.csv', '_delta_avg.csv')
outname_mdat_delta_avg = os.path.basename(mdat_f).replace(".csv", "_delta_avg.csv")

mc_all.write_csv(os.path.join(out_dir, outname_concat))
print("%s saved at %s." % (outname_concat, out_dir))
delta_all.write_csv(os.path.join(out_dir, outname_delta))
print("%s saved at %s." % (outname_delta, out_dir))
mdat_delta_all.write_csv(os.path.join(out_dir, outname_mdat_delta))
print("%s saved at %s." % (outname_mdat_delta, out_dir))

delta_avg.write_csv(os.path.join(out_dir, outname_delta_avg))
print("%s saved at %s." % (outname_delta_avg, out_dir))
mc_mdat_avg.write_csv(os.path.join(out_dir, outname_mdat_delta_avg))
print("%s saved at %s." % (outname_mdat_delta_avg, out_dir))