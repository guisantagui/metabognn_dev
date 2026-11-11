#!/usr/bin/env python
"""
Do the expression signatures carry drug-specific information? Model-free check, before any target prediction.

Active perturbations (>= --min-de-genes DE genes) of drugs with known targets are compared by the Spearman
correlation of their full signatures (all genes of --signatures, e.g. limma log2FC). If signatures say something
about the drug, perturbations of drugs that share a target should be more similar than those of unrelated drugs.

Reported (stdout and <out-dir>/summary.json; random = AUROC 0.5, enrichment 1):
  same_drug           pairs of the same drug (other dose / time) vs different drugs sharing no target: the
                      reproducibility ceiling (positive control)
  shared_target       pairs of different drugs sharing >= 1 target vs pairs sharing none, also within the same
                      time point (time is a strong confounder of signature similarity)
  shared_target_drug  the same on drug consensus signatures (mean of the drug's active perturbations)
  knn                 for each perturbation, the fraction of its --k nearest perturbations of other drugs that
                      share a target with it, over the fraction expected at random
  *_null              the same with drug labels permuted among perturbations (--n-perm times)
<out-dir>/per_target.tsv: for each target hit by >= --min-drugs drugs, AUROC of pairs where both drugs hit it
vs pairs where only one does (and they share no target): which target classes leave a recognisable signature.
<out-dir>/overlap_drug_pairs.tsv: the drug pairs behind each target-overlap bin (perturbation pairs, mean
correlation overall, at the same and at different time points, shared targets), largest first.
<out-dir>/top_pairs_different_drugs.tsv: the 300 most correlated perturbation pairs of different drugs, to spot
the same compound listed under two names.
For plot_signature_target_similarity.py: <out-dir>/pairs_sample.tsv.gz (correlations of up to --plot-pairs pairs per
class: same drug, shared target, no shared target) and example_targets_corr.tsv.gz (perturbation x perturbation
correlations of the drugs of up to --plot-targets well-recognised targets, each drug under one target only).
"""
import argparse
import json
import os

import numpy as np
import pandas as pd
from scipy import sparse

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("--signatures", required=True, help="genes x perturbations CSV (limma log2FC.csv)")
ap.add_argument("--perturbations", required=True, help="prepare_inputs.py perturbations.tsv")
ap.add_argument("--drug-targets", required=True, help="build_drug_targets.py drug_targets.tsv")
ap.add_argument("--de-counts", required=True, help="perturbation -> number of DE genes (fit_info_nsig.csv)")
ap.add_argument("--de-column", default="n_sig")
ap.add_argument("--min-de-genes", type=int, default=200)
ap.add_argument("--top-var", type=int, default=0, help="use only the N most variable genes (0: all)")
ap.add_argument("--k", type=int, default=10)
ap.add_argument("--min-drugs", type=int, default=5, help="per_target.tsv: targets of at least this many drugs")
ap.add_argument("--n-perm", type=int, default=3)
ap.add_argument("--seed", type=int, default=0)
ap.add_argument("--plot-pairs", type=int, default=50000, help="pairs_sample.tsv.gz: pairs per class")
ap.add_argument("--plot-targets", type=int, default=6, help="example_targets_corr.tsv.gz: number of targets")
ap.add_argument("--plot-perts-per-target", type=int, default=40)
ap.add_argument("--out-dir", required=True)
args = ap.parse_args()
os.makedirs(args.out_dir, exist_ok=True)
rng = np.random.default_rng(args.seed)


def auroc(score, label):
    """Mann-Whitney AUROC (ties broken arbitrarily, fine for continuous scores)."""
    label = np.asarray(label, bool)
    n1, n0 = label.sum(), (~label).sum()
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = np.empty(len(score))
    r[np.argsort(score, kind="stable")] = np.arange(1, len(score) + 1)
    return float((r[label].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


# ---- perturbations, activity, targets (drugs by name, pooled over BRD ids as in metabolic_target_prediction.py)
pt = pd.read_csv(args.perturbations, sep="\t", dtype=str, keep_default_na=False).set_index("perturbation")
brd_to_name = {b: n.lower() for n, ids in zip(pt.pertname, pt.pert_ids) for b in ids.split("|")}
dt = pd.read_csv(args.drug_targets, sep="\t", dtype=str, keep_default_na=False)
dt["name"] = dt.drug.map(brd_to_name).fillna(dt.pertname.str.lower())
targets = dt.groupby("name").target.agg(set).to_dict()
n_de = pd.read_csv(args.de_counts, index_col=0)[args.de_column]

header = pd.read_csv(args.signatures, nrows=0).columns
cols = [c for c in header[1:] if c in pt.index and n_de.get(c, 0) >= args.min_de_genes
        and pt.at[c, "pertname"].lower() in targets]
print(f"{len(cols):,} active perturbations (>= {args.min_de_genes} DE genes) of drugs with known targets")
sig = pd.read_csv(args.signatures, usecols=[header[0]] + cols, index_col=0).dropna()
if args.top_var:
    sig = sig.loc[sig.var(axis=1).nlargest(args.top_var).index]
print(f"signatures: {sig.shape[0]:,} genes")

drug = pt.loc[cols, "pertname"].str.lower().to_numpy()
time = pt.loc[cols, "timepoint"].to_numpy()
drugs, d = np.unique(drug, return_inverse=True)
genes = sorted(set().union(*(targets[x] for x in drugs)))
gpos = {g: i for i, g in enumerate(genes)}
rows = [i for i, x in enumerate(drugs) for g in targets[x]]
cols_t = [gpos[g] for x in drugs for g in targets[x]]
T = sparse.csr_matrix((np.ones(len(rows)), (rows, cols_t)), shape=(len(drugs), len(genes)))
inter = (T @ T.T).toarray()
share = inter > 0  # drugs x drugs: at least one target in common
n_t = np.asarray(T.sum(1)).ravel()
jac = inter / (n_t[:, None] + n_t[None, :] - inter)  # Jaccard of the two target sets
print(f"{len(drugs):,} drugs, {len(genes):,} distinct targets")

# ---- Spearman correlation = Pearson on per-signature gene ranks
Z = sig.rank().to_numpy(np.float32)
Z = (Z - Z.mean(0)) / Z.std(0)
C = (Z.T @ Z) / Z.shape[0]
i, j = np.triu_indices(len(cols), 1)
c = C[i, j]
same = d[i] == d[j]
shr = share[d[i], d[j]] & ~same
none = ~share[d[i], d[j]]
same_t = time[i] == time[j]


def pair_metrics(dd):
    sm = dd[i] == dd[j]
    sh = share[dd[i], dd[j]] & ~sm
    no = ~share[dd[i], dd[j]]
    m = sm | no
    o = sh | no
    ot = o & same_t
    return {"same_drug": auroc(c[m], sm[m]), "shared_target": auroc(c[o], sh[o]),
            "shared_target_same_time": auroc(c[ot], sh[ot])}


def knn_enrichment(dd):
    S = C.copy()
    S[dd[:, None] == dd[None, :]] = -np.inf  # neighbours among other drugs only
    nn = np.argpartition(-S, args.k, axis=1)[:, :args.k]
    hit = share[dd[:, None], dd[nn]].mean(1)
    other = dd[:, None] != dd[None, :]
    base = (share[dd[:, None], dd[None, :]] & other).sum(1) / other.sum(1)
    return float(hit.mean() / base.mean()), float(hit.mean()), float(base.mean())


res = {"n_perturbations": len(cols), "n_drugs": len(drugs), "n_genes": int(sig.shape[0]),
       "n_pairs": {"same_drug": int(same.sum()), "shared_target": int(shr.sum()), "no_shared_target": int(none.sum())},
       "mean_corr": {"same_drug": float(c[same].mean()), "shared_target": float(c[shr].mean()),
                     "no_shared_target": float(c[none].mean())}}
res.update(pair_metrics(d))

# graded: pairs of different drugs by the overlap of their target sets, each bin against pairs sharing none
jp = jac[d[i], d[j]]
res["by_target_overlap"] = {}
bin_rows = []
for lo, hi, name in ((0, 0.1, "jaccard (0, 0.1]"), (0.1, 0.25, "jaccard (0.1, 0.25]"), (0.25, 0.5, "jaccard (0.25, 0.5]"),
                     (0.5, 0.999, "jaccard (0.5, 1)"), (0.999, 1.0, "identical target sets")):
    pos = (jp > lo) & (jp <= hi) & ~same
    m = pos | none
    st, dt_ = pos & same_t, pos & ~same_t
    res["by_target_overlap"][name] = {"n_pairs": int(pos.sum()), "n_drug_pairs": int(len(set(zip(d[i][pos], d[j][pos])))),
                                      "mean_corr": float(c[pos].mean()) if pos.any() else None,
                                      "mean_corr_same_time": float(c[st].mean()) if st.any() else None,
                                      "mean_corr_diff_time": float(c[dt_].mean()) if dt_.any() else None,
                                      "auroc_vs_no_shared": auroc(c[m], pos[m])}
    if pos.any():
        bin_rows.append(pd.DataFrame({"overlap": name, "a": d[i][pos], "b": d[j][pos], "corr": c[pos],
                                      "same_time": same_t[pos]}))

# drug pairs behind each overlap bin: which ones dominate it (many perturbation pairs) and whether time drives them
bp = pd.concat(bin_rows)
lo_, hi_ = np.minimum(bp.a, bp.b), np.maximum(bp.a, bp.b)
bp["drug_a"], bp["drug_b"] = drugs[lo_], drugs[hi_]
g = bp.groupby(["overlap", "drug_a", "drug_b"], sort=False)
drug_pairs = g.agg(n_pert_pairs=("corr", "size"), mean_corr=("corr", "mean"),
                   frac_same_time=("same_time", "mean")).reset_index()
drug_pairs["mean_corr_same_time"] = g.apply(lambda x: x["corr"][x.same_time].mean()).to_numpy()
drug_pairs["mean_corr_diff_time"] = g.apply(lambda x: x["corr"][~x.same_time].mean()).to_numpy()
drug_pairs["shared_targets"] = [",".join(sorted(targets[a] & targets[b])[:10]) for a, b in
                                zip(drug_pairs.drug_a, drug_pairs.drug_b)]
drug_pairs.sort_values(["overlap", "n_pert_pairs"], ascending=[True, False]).to_csv(
    os.path.join(args.out_dir, "overlap_drug_pairs.tsv"), sep="\t", index=False, float_format="%.4g")

# most correlated perturbation pairs of different drugs: duplicated compounds under two names show up here
dose = pt.loc[cols, "dose"].to_numpy()
other = ~same
top = np.argsort(-np.where(other, c, -np.inf))[:300]
pd.DataFrame({"pert_a": np.array(cols)[i[top]], "pert_b": np.array(cols)[j[top]], "drug_a": drug[i[top]],
              "drug_b": drug[j[top]], "dose_a": dose[i[top]], "dose_b": dose[j[top]], "time_a": time[i[top]],
              "time_b": time[j[top]], "corr": c[top], "target_jaccard": jp[top],
              "identical_targets": jp[top] > 0.999}).to_csv(
    os.path.join(args.out_dir, "top_pairs_different_drugs.tsv"), sep="\t", index=False, float_format="%.4g")
enr, hit, base = knn_enrichment(d)
res["knn"] = {"k": args.k, "enrichment": enr, "frac_sharing_target": hit, "expected": base}

# drug consensus signatures
D = np.stack([Z[:, d == q].mean(1) for q in range(len(drugs))], 1)
D = (D - D.mean(0)) / D.std(0)
CD = (D.T @ D) / D.shape[0]
a, b = np.triu_indices(len(drugs), 1)
res["shared_target_drug"] = auroc(CD[a, b], share[a, b])

# null: drug labels permuted among perturbations (each permutation keeps the number of perturbations per drug)
nulls = [pair_metrics(rng.permutation(d)) | {"knn_enrichment": knn_enrichment(rng.permutation(d))[0]}
         for _ in range(args.n_perm)]
res["null_mean"] = {k: float(np.nanmean([n[k] for n in nulls])) for k in nulls[0]}

# ---- per target
tr = []
for g, col in gpos.items():
    has = T[:, col].toarray().ravel() > 0
    if has.sum() < args.min_drugs:
        continue
    hi, hj = has[d[i]], has[d[j]]
    pos = hi & hj & ~same
    neg = (hi ^ hj) & none
    m = pos | neg
    if pos.sum() and neg.sum():
        tr.append(dict(target=g, n_drugs=int(has.sum()), n_perts=int(has[d].sum()), n_pos_pairs=int(pos.sum()),
                       auroc=auroc(c[m], pos[m]), mean_corr_pos=float(c[pos].mean()),
                       mean_corr_neg=float(c[neg].mean())))
per_target = pd.DataFrame(tr).sort_values("auroc", ascending=False)
per_target.to_csv(os.path.join(args.out_dir, "per_target.tsv"), sep="\t", index=False)
with open(os.path.join(args.out_dir, "summary.json"), "w") as fh:
    json.dump(res, fh, indent=2)

# ---- for plot_signature_target_similarity.py
samp = []
for name, m in (("same drug", same), ("shared target", shr), ("no shared target", none)):
    idx = np.flatnonzero(m)
    if len(idx) > args.plot_pairs:
        idx = np.sort(rng.choice(idx, args.plot_pairs, replace=False))
    samp.append(pd.DataFrame({"pair_class": name, "corr": c[idx], "target_jaccard": jp[idx], "same_time": same_t[idx]}))
pd.concat(samp).to_csv(os.path.join(args.out_dir, "pairs_sample.tsv.gz"), sep="\t", index=False, float_format="%.4g")

# example block: the best-recognised targets (per_target order), each drug shown under the first target it falls in
ex, used = [], set()
for t in per_target.target if len(per_target) else []:
    if len(ex) == args.plot_targets:
        break
    tdrugs = [q for q in np.flatnonzero(T[:, gpos[t]].toarray().ravel()) if q not in used]
    if len(tdrugs) < args.min_drugs:
        continue
    used.update(tdrugs)
    p = np.flatnonzero(np.isin(d, tdrugs))
    if len(p) > args.plot_perts_per_target:
        p = np.sort(rng.choice(p, args.plot_perts_per_target, replace=False))
    hours = pd.to_numeric(pd.Series(time[p]).str.extract(r"([\d.]+)")[0], errors="coerce").to_numpy()
    ex.append(pd.DataFrame({"target": t, "drug": drug[p], "timepoint": time[p], "perturbation": np.array(cols)[p],
                            "_h": hours, "_k": p}).sort_values(["drug", "_h"]))
if ex:
    ex = pd.concat(ex, ignore_index=True)
    block = pd.DataFrame(C[np.ix_(ex._k, ex._k)], columns=ex.perturbation)
    pd.concat([ex.drop(columns=["_h", "_k"]), block], axis=1).to_csv(
        os.path.join(args.out_dir, "example_targets_corr.tsv.gz"), sep="\t", index=False, float_format="%.3g")

print(json.dumps(res, indent=2))
if len(per_target):
    print(f"\nper target (>= {args.min_drugs} drugs), top and bottom by AUROC:")
    print(pd.concat([per_target.head(15), per_target.tail(5)]).to_string(index=False, float_format="%.3f"))
print(f"\nsaved {args.out_dir}/summary.json, per_target.tsv, pairs_sample.tsv.gz, example_targets_corr.tsv.gz; "
      f"plots: python plot_signature_target_similarity.py --dir {args.out_dir}")
