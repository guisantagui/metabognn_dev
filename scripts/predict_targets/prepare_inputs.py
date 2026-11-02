#!/usr/bin/env python
"""
Inputs of metabolic_target_prediction.py that only need computing once per cell line.

  baseline_ctrl.tsv   symbol, value: mean expression over the control samples (exper_group == "ctrl"), with LINCS
                      symbols mapped to HGNC approved symbols (aliases and previous symbols resolved; genes mapping
                      to the same approved symbol are averaged; genes without an HGNC protein-coding symbol kept
                      under their LINCS name, they simply match no network node)
  perturbations.tsv   perturbation, pertname, pert_ids, dose, timepoint, n_samples (treated samples only)

Usage
  python prepare_inputs.py --cell-dir ../../data/perturbation/lincs/cell_matrices/A549 \
                           --hgnc ../../data/utilities/hgnc_complete_set.txt
"""
import argparse
import os
import sys
import time

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "build_network"))
from build_network import HGNC  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cell-dir", required=True, help="folder with expr.csv (genes x samples) and mdat.csv")
    ap.add_argument("--hgnc", required=True)
    ap.add_argument("--out-dir", help="defaults to --cell-dir")
    args = ap.parse_args()
    out_dir = args.out_dir or args.cell_dir

    mdat = pd.read_csv(os.path.join(args.cell_dir, "mdat.csv"), index_col=0, dtype=str, keep_default_na=False)
    ctrl = mdat.id[mdat.exper_group == "ctrl"].tolist()
    t0 = time.time()
    f = os.path.join(args.cell_dir, "expr.csv")
    try:  # polars reads only the requested columns quickly; pandas parses every line in full
        import polars as pl
        expr = pl.read_csv(f, columns=["V1"] + ctrl).to_pandas().set_index("V1")
    except ImportError:
        expr = pd.read_csv(f, usecols=["V1"] + ctrl, index_col="V1")
    print(f"{expr.shape[0]:,} genes x {expr.shape[1]:,} control samples read in {time.time() - t0:.0f}s")

    hgnc = HGNC(args.hgnc)
    approved = expr.index.map(lambda s: hgnc(s) or s)
    print(f"{(approved != expr.index).sum():,} LINCS symbols renamed to the HGNC approved symbol")
    base = expr.mean(axis=1).groupby(approved).mean()
    base.rename_axis("symbol").rename("value").to_csv(os.path.join(out_dir, "baseline_ctrl.tsv"), sep="\t")
    print(f"baseline_ctrl.tsv: {len(base):,} genes")

    pert = mdat[mdat.exper_group != "ctrl"]
    table = pert.groupby("perturbation").agg(pertname=("pertname", "first"),
                                             pert_ids=("pert_id", lambda s: "|".join(sorted(set(s)))),
                                             dose=("dose", "first"), timepoint=("timepoint", "first"),
                                             n_samples=("id", "size")).reset_index()
    table.to_csv(os.path.join(out_dir, "perturbations.tsv"), sep="\t", index=False)
    print(f"perturbations.tsv: {len(table):,} perturbations of {table.pertname.nunique():,} drugs")


if __name__ == "__main__":
    main()
