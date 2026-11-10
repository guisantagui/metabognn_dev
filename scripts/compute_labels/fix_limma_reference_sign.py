#!/usr/bin/env python
"""
Correct limma_DE.R outputs written before the control reference level was set explicitly.

With design = model.matrix(~ perturbation), the reference level was whichever of the perturbation name and the
control name ("DMSO_nan_<time>") sorted first in the R session's locale. Where the perturbation came first,
limma returned log2FC = control - treated and basal = the treated mean. p- and q-values are unaffected
(two-sided), so the fix is: log2FC -> -log2FC and basal -> basal + log2FC (the control mean).

Which perturbations were affected is decided from the data, not from the name (the collation of the locale the job
ran in is unknown here): for each perturbation the two candidates for its control mean, basal (not flipped) and
basal + log2FC (flipped), are compared with a reference control profile (per time point, the feature-wise median
of the current control estimates, refined twice). Perturbations whose two candidates are about as close to the
reference (no effect, so no way to tell) follow the name rule (case-insensitive alphanumeric comparison with the
control name, as en_US-type collations do); its agreement with the data on the clear cases is reported, and if it
is not ~100% the decision should not be trusted.

Writes <out-dir>/log2FC.csv, basal.csv and flipped.tsv (per perturbation: decision, evidence, name rule).
qvals.csv, pvals.csv and fit_info.csv do not change: copy or link them.
"""
import argparse
import os

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--dir", required=True, help="limma_DE.R output folder (log2FC.csv, basal.csv)")
ap.add_argument("--out-dir", required=True)
ap.add_argument("--control", default="DMSO_nan", help="control name without the time point")
ap.add_argument("--min-ratio", type=float, default=1.5,
                help="clear case: the farther candidate is at least this many times farther from the reference")
args = ap.parse_args()
if os.path.abspath(args.out_dir) == os.path.abspath(args.dir):
    ap.error("--out-dir must differ from --dir (the originals are kept)")
os.makedirs(args.out_dir, exist_ok=True)

fc = pd.read_csv(os.path.join(args.dir, "log2FC.csv"), index_col=0, dtype={0: str}).astype(np.float32)
basal = pd.read_csv(os.path.join(args.dir, "basal.csv"), index_col=0, dtype={0: str}).astype(np.float32)
if not (fc.index.equals(basal.index) and fc.columns.equals(basal.columns)):
    raise SystemExit("log2FC.csv and basal.csv have different rows or columns")
print(f"{fc.shape[0]:,} features x {fc.shape[1]:,} perturbations")

time = fc.columns.to_series().str.extract(r"_([^_]*)$")[0].fillna("")
A = basal.to_numpy()                       # control mean if not flipped
B = A + np.nan_to_num(fc.to_numpy())       # control mean if flipped
flip = np.zeros(fc.shape[1], bool)
dA = dB = None
for _ in range(3):
    ctrl = np.where(flip, B, A)
    dA, dB = np.empty(len(flip)), np.empty(len(flip))
    for t in time.unique():
        cols = np.where(time.to_numpy() == t)[0]
        ref = np.nanmedian(ctrl[:, cols], axis=1)
        dA[cols] = np.sqrt(np.nanmean((A[:, cols] - ref[:, None]) ** 2, axis=0))
        dB[cols] = np.sqrt(np.nanmean((B[:, cols] - ref[:, None]) ** 2, axis=0))
    flip = dB < dA


def key(s):
    return "".join(ch for ch in s.lower() if ch.isalnum())


name_rule = np.array([key(p) < key(f"{args.control}_{t}") for p, t in zip(fc.columns, time)])
ratio = np.maximum(dA, dB) / np.maximum(np.minimum(dA, dB), 1e-12)
clear = ratio >= args.min_ratio
agree = (flip == name_rule)[clear].mean() if clear.any() else float("nan")
final = np.where(clear, flip, name_rule)
print(f"clear cases: {clear.sum():,} of {len(flip):,} (distance ratio >= {args.min_ratio}); "
      f"flipped among them {flip[clear].sum():,}; name rule agrees on {agree:.2%}")
print(f"undecidable from the data (no effect): {(~clear).sum():,}, set by the name rule "
      f"({name_rule[~clear].sum():,} flipped)")
print(f"total flipped: {final.sum():,} ({final.mean():.1%})")
if agree < 0.98:
    print("WARNING: the name rule disagrees with the data on clear cases; check flipped.tsv before using the output")

pd.DataFrame({"perturbation": fc.columns, "flipped": final, "flipped_by_data": flip, "clear": clear,
              "distance_ratio": ratio, "name_rule": name_rule}).to_csv(
    os.path.join(args.out_dir, "flipped.tsv"), sep="\t", index=False, float_format="%.4g")
fc.loc[:, final] = -fc.loc[:, final]
basal.loc[:, final] = B[:, final]
fc.to_csv(os.path.join(args.out_dir, "log2FC.csv"), float_format="%.6g")
basal.to_csv(os.path.join(args.out_dir, "basal.csv"), float_format="%.6g")
print(f"saved {args.out_dir}/log2FC.csv, basal.csv, flipped.tsv")
