#!/usr/bin/env python
"""
Check fix_limma_reference_sign.py against limma itself, on a few perturbations.

  select   pick perturbations from <fixed-dir>/flipped.tsv: the clear cases where the data disagrees with R's
           level order (the doubtful ones, up to --n-disagree), plus --n-control flipped and --n-control not flipped
           where both agree; writes a list for limma_DE.R --perts
  compare  compare the corrected log2FC.csv / basal.csv with a limma_DE.R --perts rerun (fixed limma_DE.R, explicit
           control reference): per perturbation the correlation and largest difference of log2FC, and which of
           the two decisions (R order, data) the rerun confirms

Typical use:
  python check_limma_sign_fix.py select --fixed-dir <expr_signfixed> --out perts_check.txt
  Rscript limma_DE.R --expr <expr_ensembl.csv> --mdat <mdat.csv> --workers 8 --out_dir <check_dir> --perts perts_check.txt
  python check_limma_sign_fix.py compare --fixed-dir <expr_signfixed> --rerun-dir <check_dir>
"""
import argparse
import os

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("mode", choices=["select", "compare"])
ap.add_argument("--fixed-dir", required=True, help="fix_limma_reference_sign.py output (flipped.tsv, log2FC.csv, basal.csv)")
ap.add_argument("--rerun-dir", help="compare: limma_DE.R --perts output")
ap.add_argument("--out", default="perts_check.txt", help="select: list of perturbations")
ap.add_argument("--n-disagree", type=int, default=30)
ap.add_argument("--n-control", type=int, default=10)
ap.add_argument("--seed", type=int, default=0)
args = ap.parse_args()

fl = pd.read_csv(os.path.join(args.fixed_dir, "flipped.tsv"), sep="\t").set_index("perturbation")
for c in ("flipped", "flipped_by_data", "clear", "name_rule"):
    fl[c] = fl[c].astype(str).str.upper() == "TRUE"

if args.mode == "select":
    rng = np.random.default_rng(args.seed)
    dis = fl[fl.clear & (fl.flipped_by_data != fl.name_rule)].sort_values("distance_ratio", ascending=False)
    ok = fl[fl.clear & (fl.flipped_by_data == fl.name_rule)]
    pick = list(dis.index[:args.n_disagree])
    for f in (True, False):
        pool = ok.index[ok.name_rule == f]
        pick += list(rng.choice(pool, min(args.n_control, len(pool)), replace=False))
    with open(args.out, "w") as fh:
        fh.write("\n".join(pick) + "\n")
    print(f"{len(pick)} perturbations ({min(len(dis), args.n_disagree)} where data and R order disagree, "
          f"{len(pick) - min(len(dis), args.n_disagree)} where they agree) written to {args.out}")
    raise SystemExit

if not args.rerun_dir:
    ap.error("compare needs --rerun-dir")
ref = pd.read_csv(os.path.join(args.rerun_dir, "log2FC.csv"), index_col=0)
perts = list(ref.columns)
f_fixed = os.path.join(args.fixed_dir, "log2FC.csv")
idx_col = pd.read_csv(f_fixed, nrows=0).columns[0]
fixed = pd.read_csv(f_fixed, index_col=0, usecols=[idx_col] + perts).loc[ref.index, perts]

rows = []
for p in perts:
    a, b = ref[p].to_numpy(), fixed[p].to_numpy()
    ok = ~(np.isnan(a) | np.isnan(b))
    r = np.corrcoef(a[ok], b[ok])[0, 1]
    # the rerun is correct by construction; r ~ +1: the applied decision was right, r ~ -1: it was wrong
    applied_right = r > 0
    rows.append(dict(perturbation=p, corr_with_rerun=r, max_abs_diff=float(np.max(np.abs(a[ok] - b[ok]))),
                     applied_flip=fl.at[p, "flipped"], r_order=fl.at[p, "name_rule"], data=fl.at[p, "flipped_by_data"],
                     clear=fl.at[p, "clear"], distance_ratio=fl.at[p, "distance_ratio"],
                     truth_flipped=fl.at[p, "flipped"] if applied_right else not fl.at[p, "flipped"]))
res = pd.DataFrame(rows)
res["r_order_right"] = res.r_order == res.truth_flipped
res["data_right"] = res.data == res.truth_flipped
res.to_csv(os.path.join(args.rerun_dir, "sign_fix_check.tsv"), sep="\t", index=False, float_format="%.4g")
dis = res.r_order != res.data
print(f"{len(res)} perturbations; corrected log2FC vs rerun: correlation median {res.corr_with_rerun.median():.4f}, "
      f"{(res.corr_with_rerun > 0.99).mean():.0%} above 0.99, {(res.corr_with_rerun < 0).sum()} negative")
print(f"where R order and data disagree ({dis.sum()}): R order right {res.r_order_right[dis].mean():.0%}, "
      f"data right {res.data_right[dis].mean():.0%}")
print(f"where they agree ({(~dis).sum()}): right {res.r_order_right[~dis].mean():.0%}")
print(f"saved {os.path.join(args.rerun_dir, 'sign_fix_check.tsv')}")
