#!/usr/bin/env python
"""
Correct limma_DE.R outputs written before the control reference level was set explicitly.

With design = model.matrix(~ perturbation), the reference level was whichever of the perturbation name and the
control name ("DMSO_nan_<time>") sorted first in the R session's locale. Where the perturbation came first,
limma returned log2FC = control - treated and basal = the treated mean. p- and q-values are unaffected
(two-sided), so the fix is: log2FC -> -log2FC and basal -> basal + log2FC (the control mean).

With --name-flips (limma_reference_order.R run in the limma jobs' environment) R's own level order decides every
perturbation, which is exactly what limma did, and the data below is only a check. Without it, which
perturbations were affected is decided from the data, not from the name (the collation of the locale the job
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
ap.add_argument("--force", action="store_true",
                help="apply --name-flips even when the data says the outputs are not flipped")
ap.add_argument("--block", type=int, default=500, help="features per block when deciding the flips")
ap.add_argument("--name-flips", help="limma_reference_order.R output (perturbation, flipped_R): R's own level order "
                                     "in the limma jobs' environment, used instead of the approximate name rule")
args = ap.parse_args()
if os.path.abspath(args.out_dir) == os.path.abspath(args.dir):
    ap.error("--out-dir must differ from --dir (the originals are kept)")
if os.path.exists(os.path.join(args.dir, "reference_level.txt")):
    raise SystemExit(f"{args.dir} was written by the fixed limma_DE.R (reference_level.txt): nothing to correct")
os.makedirs(args.out_dir, exist_ok=True)

def read(name):
    """float32 straight from the parser (no float64 copy): memory ~ 4 bytes per value."""
    f = os.path.join(args.dir, name)
    cols = pd.read_csv(f, nrows=0).columns
    df = pd.read_csv(f, index_col=0, dtype={cols[0]: str, **{c: np.float32 for c in cols[1:]}})
    return df.index, df.columns, np.ascontiguousarray(df.to_numpy(np.float32))


rows, cols, F = read("log2FC.csv")
rows_b, cols_b, Bas = read("basal.csv")
if not (rows.equals(rows_b) and cols.equals(cols_b)):
    raise SystemExit("log2FC.csv and basal.csv have different rows or columns")
print(f"{F.shape[0]:,} features x {F.shape[1]:,} perturbations")

# Candidates for each perturbation's control mean: A = basal (not flipped), B = basal + log2FC (flipped). They are
# computed in blocks of features, so that no other full-size matrix is held besides the two inputs.
time = cols.to_series().str.extract(r"_([^_]*)$")[0].fillna("").to_numpy()
groups = [np.where(time == t)[0] for t in np.unique(time)]
flip = np.zeros(F.shape[1], bool)
dA = dB = None
for _ in range(3):
    sA, sB, n = np.zeros(F.shape[1]), np.zeros(F.shape[1]), np.zeros(F.shape[1])
    for r in range(0, F.shape[0], args.block):
        A = Bas[r:r + args.block].astype(np.float64)
        B = A + np.nan_to_num(F[r:r + args.block])
        ctrl = np.where(flip, B, A)
        for g in groups:
            ref = np.nanmedian(ctrl[:, g], axis=1)[:, None]
            ea, eb = (A[:, g] - ref) ** 2, (B[:, g] - ref) ** 2
            sA[g] += np.nansum(ea, 0)
            sB[g] += np.nansum(eb, 0)
            n[g] += (~np.isnan(ea)).sum(0)
    dA, dB = np.sqrt(sA / np.maximum(n, 1)), np.sqrt(sB / np.maximum(n, 1))
    flip = dB < dA


def key(s):
    return "".join(ch for ch in s.lower() if ch.isalnum())


if args.name_flips:
    fr = pd.read_csv(args.name_flips, sep="\t").set_index("perturbation").flipped_R.astype(str).str.upper() == "TRUE"
    missing = ~cols.isin(fr.index)
    if missing.any():
        raise SystemExit(f"{missing.sum():,} perturbations not in --name-flips, e.g. {list(cols[missing][:3])}")
    name_rule = fr.reindex(cols).to_numpy()
    print(f"name rule: R's level order from {args.name_flips}")
else:
    name_rule = np.array([key(p) < key(f"{args.control}_{t}") for p, t in zip(cols, time)])
ratio = np.maximum(dA, dB) / np.maximum(np.minimum(dA, dB), 1e-12)
clear = ratio >= args.min_ratio
agree = (flip == name_rule)[clear].mean() if clear.any() else float("nan")
print(f"clear cases: {clear.sum():,} of {len(flip):,} (distance ratio >= {args.min_ratio}); "
      f"flipped among them {flip[clear].sum():,}; name rule agrees on {agree:.2%}")
if args.name_flips and clear.sum() >= 100 and flip[clear].mean() < 0.05 and name_rule.mean() > 0.2 and not args.force:
    # e.g. outputs of the fixed limma_DE.R written before the marker existed: applying R's (old) order would flip
    # correct signatures
    raise SystemExit(f"the data finds almost no flipped perturbation ({flip[clear].mean():.1%} of the clear cases) "
                     f"while R's level order would flip {name_rule.mean():.1%}: these outputs look already correct "
                     f"(fixed limma_DE.R). Nothing written; --force to apply anyway")
if args.name_flips:
    # R's level order is what limma did: it decides every perturbation, the data is only a check
    final = name_rule
    print(f"decision: R's level order for all perturbations; the data disagrees on {(flip != name_rule)[clear].sum():,} "
          f"clear cases (listed in flipped.tsv: clear & flipped_by_data != flipped)")
else:
    final = np.where(clear, flip, name_rule)
    print(f"undecidable from the data (no effect): {(~clear).sum():,}, set by the name rule "
          f"({name_rule[~clear].sum():,} flipped)")
    if agree < 0.98:
        print("WARNING: the name rule disagrees with the data on clear cases; check flipped.tsv before using the output")
print(f"total flipped: {final.sum():,} ({final.mean():.1%})")

pd.DataFrame({"perturbation": cols, "flipped": final, "flipped_by_data": flip, "clear": clear,
              "distance_ratio": ratio, "name_rule": name_rule}).to_csv(
    os.path.join(args.out_dir, "flipped.tsv"), sep="\t", index=False, float_format="%.4g")
# in place, one flipped column at a time: basal <- basal + log2FC (the control mean), log2FC <- -log2FC
for k in np.where(final)[0]:
    Bas[:, k] += np.nan_to_num(F[:, k])
    F[:, k] *= -1
for name, X in (("log2FC.csv", F), ("basal.csv", Bas)):
    pd.DataFrame(X, index=rows, columns=cols, copy=False).to_csv(os.path.join(args.out_dir, name),
                                                                 float_format="%.6g", chunksize=1000)
print(f"saved {args.out_dir}/log2FC.csv, basal.csv, flipped.tsv")
