#!/usr/bin/env python
"""
Compare two complete limma label directories (log2FC.csv, qvals.csv, optionally basal.csv / pvals.csv), e.g. the
copy corrected by fix_limma_reference_sign.py against a rerun of the fixed limma_DE.R.

Files are features x perturbations (first column = feature id). Rows and columns are matched by name, not position.
Each matrix is read in row chunks as float32 (expr A549: ~20k x ~6k -> ~0.5 GB per matrix; one matrix pair in
memory at a time).

Per perturbation (written to <out>/per_perturbation.tsv):
  log2FC  Pearson r, max |A - B|, fraction of entries with opposite sign  -> status identical / negated / differs
  qvals   max |A - B|, number with q < --q in A and B, Jaccard of those sets (this drives the activity filter)
  top-n   overlap of the top --top-n features by |log2FC| and sign agreement on the shared ones: what the model reads
          under --ignore-qvalue --top-n (approximation: the model ranks only features mapped to network nodes)
  basal / pvals  max |A - B| (if present in both)

Run on the cluster, where both directories are on local disk (reading these files over a network mount can hang).
  python compare_label_dirs.py --a <dir_corrected> --b <dir_rerun> --out <out_dir> [--max-perts 200]
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--a", required=True, help="label directory A (e.g. sign-corrected)")
ap.add_argument("--b", required=True, help="label directory B (e.g. rerun with fixed limma_DE.R)")
ap.add_argument("--out", required=True, help="output directory")
ap.add_argument("--tol", type=float, default=1e-4, help="max |A - B| still counted as identical (default 1e-4)")
ap.add_argument("--q", type=float, default=0.05, help="q-value threshold for significance counts")
ap.add_argument("--top-n", type=int, default=100)
ap.add_argument("--max-perts", type=int, default=0, help="compare a random subset of common perturbations (0 = all)")
ap.add_argument("--chunksize", type=int, default=2000, help="rows per read chunk")
ap.add_argument("--seed", type=int, default=0)
args = ap.parse_args()

T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


def header(path):
    cols = pd.read_csv(path, nrows=0).columns
    return cols[0], list(cols[1:])


def load(path, cols):
    """Read the given perturbation columns of a features x perturbations CSV as float32, in row chunks."""
    idx_col = pd.read_csv(path, nrows=0).columns[0]
    parts, n = [], 0
    for ch in pd.read_csv(path, index_col=0, usecols=[idx_col] + cols, chunksize=args.chunksize,
                          dtype={c: np.float32 for c in cols}):
        parts.append(ch)
        n += len(ch)
        if len(parts) % 5 == 0:
            log(f"  {os.path.basename(path)}: {n:,} rows")
    df = pd.concat(parts)
    df.index = df.index.astype(str)
    return df


def paths(name):
    return os.path.join(args.a, name), os.path.join(args.b, name)


def align(name, perts):
    pa, pb = paths(name)
    log(f"reading {name} from A")
    A = load(pa, perts)
    log(f"reading {name} from B")
    B = load(pb, perts)
    dup_a, dup_b = A.index.duplicated().sum(), B.index.duplicated().sum()
    if dup_a or dup_b:
        log(f"  WARNING {name}: duplicated feature ids (A {dup_a}, B {dup_b}); first occurrence kept")
        A, B = A[~A.index.duplicated()], B[~B.index.duplicated()]
    only_a, only_b = A.index.difference(B.index), B.index.difference(A.index)
    common = A.index.intersection(B.index)
    log(f"  {name}: {len(A):,} / {len(B):,} features, {len(common):,} common, {len(only_a):,} only in A, "
        f"{len(only_b):,} only in B" + (f" (e.g. A: {list(only_a[:3])}, B: {list(only_b[:3])})" if len(only_a) or len(only_b) else ""))
    same_order = len(only_a) == 0 and len(only_b) == 0 and (A.index == B.index).all()
    if not same_order and not (len(only_a) or len(only_b)):
        log(f"  {name}: same features, different row order (aligned by name)")
    return A.loc[common, perts].to_numpy(), B.loc[common, perts].to_numpy(), common


os.makedirs(args.out, exist_ok=True)
for d in (args.a, args.b):
    for f in ("log2FC.csv", "qvals.csv"):
        if not os.path.isfile(os.path.join(d, f)):
            sys.exit(f"missing {os.path.join(d, f)}")

# perturbations
_, pa = header(os.path.join(args.a, "log2FC.csv"))
_, pb = header(os.path.join(args.b, "log2FC.csv"))
perts = [p for p in pa if p in set(pb)]
only_a, only_b = sorted(set(pa) - set(pb)), sorted(set(pb) - set(pa))
log(f"perturbations: A {len(pa):,}, B {len(pb):,}, common {len(perts):,}, only A {len(only_a):,}, only B {len(only_b):,}")
for p in ("DMSO",):
    na, nb = sum(p in c for c in only_a), sum(p in c for c in only_b)
    if na or nb:
        log(f"  of these, containing '{p}': only A {na}, only B {nb}")
if only_a or only_b:
    pd.Series(only_a, dtype=str).to_csv(os.path.join(args.out, "perts_only_in_A.txt"), index=False, header=False)
    pd.Series(only_b, dtype=str).to_csv(os.path.join(args.out, "perts_only_in_B.txt"), index=False, header=False)
if args.max_perts and len(perts) > args.max_perts:
    rng = np.random.default_rng(args.seed)
    perts = list(rng.choice(perts, args.max_perts, replace=False))
    log(f"random subset of {len(perts)} perturbations (--max-perts)")

res = pd.DataFrame(index=pd.Index(perts, name="perturbation"))

# log2FC
A, B, feats = align("log2FC.csv", perts)
nan_a, nan_b = np.isnan(A), np.isnan(B)
res["lfc_nan_mismatch"] = (nan_a != nan_b).sum(0)
ok = ~(nan_a | nan_b)
r, mad, opp = [], [], []
for j in range(len(perts)):
    a, b = A[ok[:, j], j].astype(np.float64), B[ok[:, j], j].astype(np.float64)
    r.append(np.corrcoef(a, b)[0, 1] if len(a) > 2 and a.std() > 0 and b.std() > 0 else np.nan)
    mad.append(np.abs(a - b).max() if len(a) else np.nan)
    nz = (np.abs(a) > 1e-8) & (np.abs(b) > 1e-8)
    opp.append(np.mean(np.sign(a[nz]) != np.sign(b[nz])) if nz.any() else np.nan)
res["lfc_r"], res["lfc_max_abs_diff"], res["lfc_frac_opposite_sign"] = r, mad, opp

# top-n by |log2FC| (rank method "first" in the model; compared as sets here)
ov, sgn = [], []
for j in range(len(perts)):
    a = np.where(nan_a[:, j], -np.inf, np.abs(A[:, j]))
    b = np.where(nan_b[:, j], -np.inf, np.abs(B[:, j]))
    ta = np.argsort(-a, kind="stable")[:args.top_n]
    tb = np.argsort(-b, kind="stable")[:args.top_n]
    shared = np.intersect1d(ta, tb)
    ov.append(len(shared) / args.top_n)
    sgn.append(np.mean(np.sign(A[shared, j]) == np.sign(B[shared, j])) if len(shared) else np.nan)
res[f"top{args.top_n}_overlap"], res[f"top{args.top_n}_sign_agree"] = ov, sgn
del A, B, nan_a, nan_b, ok

# qvals
A, B, _ = align("qvals.csv", perts)
with np.errstate(invalid="ignore"):
    d = np.abs(A - B)
    res["q_max_abs_diff"] = np.nanmax(np.where(np.isnan(d), -np.inf, d), axis=0)
    sa, sb = A < args.q, B < args.q
res["n_sig_A"], res["n_sig_B"] = sa.sum(0), sb.sum(0)
union = (sa | sb).sum(0)
res["sig_jaccard"] = np.where(union > 0, (sa & sb).sum(0) / np.maximum(union, 1), 1.0)
del A, B, d, sa, sb

# optional files
for name in ("basal.csv", "pvals.csv"):
    fa, fb = paths(name)
    if not (os.path.isfile(fa) and os.path.isfile(fb)):
        log(f"{name}: not in both directories, skipped")
        continue
    A, B, _ = align(name, perts)
    with np.errstate(invalid="ignore"):
        d = np.abs(A - B)
    res[f"{name[:-4]}_max_abs_diff"] = np.nanmax(np.where(np.isnan(d), -np.inf, d), axis=0)
    del A, B, d

# status per perturbation
res["status"] = "differs"
res.loc[res.lfc_max_abs_diff <= args.tol, "status"] = "identical"
res.loc[(res.status != "identical") & (res.lfc_r < -0.99), "status"] = "negated"
res.loc[(res.status == "differs") & (res.lfc_r > 0.999), "status"] = "near_identical"
out = os.path.join(args.out, "per_perturbation.tsv")
res.to_csv(out, sep="\t", float_format="%.6g")


def q(s, qs=(0, 0.5, 0.99, 1)):
    s = s.dropna()
    return " / ".join(f"{s.quantile(x):.4g}" for x in qs) if len(s) else "n/a"


print()
print(f"=== {len(res):,} perturbations compared (A = {args.a}, B = {args.b}) ===")
print(res.status.value_counts().to_string())
print(f"(identical: max |dlog2FC| <= {args.tol}; near_identical: r > 0.999; negated: r < -0.99)")
print("                           min / median / 99% / max")
print(f"log2FC r                   {q(res.lfc_r)}")
print(f"log2FC max |A-B|           {q(res.lfc_max_abs_diff)}")
print(f"log2FC frac opposite sign  {q(res.lfc_frac_opposite_sign)}")
print(f"log2FC NaN mismatches      {q(res.lfc_nan_mismatch.astype(float))}")
print(f"qvals max |A-B|            {q(res.q_max_abs_diff)}")
print(f"sig set Jaccard (q<{args.q})   {q(res.sig_jaccard)}")
print(f"n_sig A - B                {q((res.n_sig_A - res.n_sig_B).astype(float))}")
print(f"top{args.top_n} overlap              {q(res[f'top{args.top_n}_overlap'])}")
print(f"top{args.top_n} sign agreement       {q(res[f'top{args.top_n}_sign_agree'])}")
for c in ("basal_max_abs_diff", "pvals_max_abs_diff"):
    if c in res:
        print(f"{c:<27}{q(res[c])}")
bad = res[res.status != "identical"].sort_values("lfc_r")
if len(bad):
    print(f"\nworst {min(10, len(bad))} non-identical perturbations:")
    print(bad.head(10)[["lfc_r", "lfc_max_abs_diff", "lfc_frac_opposite_sign", "sig_jaccard",
                        f"top{args.top_n}_overlap", f"top{args.top_n}_sign_agree"]].to_string(float_format="%.4g"))
print(f"\nsaved {out}")
log("done")
