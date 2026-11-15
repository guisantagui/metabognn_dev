#!/usr/bin/env python
"""
Does a signature's specificity depend on how much of the perturbation's expression change falls on metabolic genes?

Perturbations are called active from transcriptome-wide DE counts (--min-de-genes), but MetChange only sees the genes
of Human-GEM (its GPR genes). A perturbation can be active with few DE genes in the GEM, leaving its MetChange signature
near noise. This script counts, per perturbation, DE genes (q < --q) overall and among the Human-GEM genes, and splits
the leak-free spec_<model>_timedose_heldout of each run and its shuffled control by that coverage. Per bin: drugs as
units (mean over the drug's perturbations in the bin), paired real − control gap, Wilcoxon p.

  python check_metabolic_coverage.py --qvals <expr labels>/qvals.csv --gem-genes Human-GEM/model/genes.tsv \
      --run mc=<real out_dir> --run expr=<real out_dir> [--out-dir DIR]

Use the same expression labels the runs used (on the cluster: the scratch copy). Each --run is a real run; its
control is <out_dir>_shuff. The spec percentiles are computed with compare_runs.py's own
functions (same strata, mid-rank percentiles), so they match compare_real_vs_shuffled/. Writes coverage.tsv (per
perturbation counts), spec_by_coverage.tsv (per run, binning and bin) and spearman.tsv.
"""
import argparse
import ast
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))


def compare_runs_functions():
    """design() and stratified_spec() from compare_runs.py (that script runs its comparison at import)."""
    src = open(os.path.join(HERE, "compare_runs.py")).read()
    tree = ast.parse(src)
    keep = [n for n in tree.body if isinstance(n, (ast.Import, ast.ImportFrom))
            or (isinstance(n, ast.FunctionDef) and n.name in ("design", "stratified_spec"))
            or (isinstance(n, ast.Assign) and any(getattr(t, "id", "") in ("STRATA", "PCT_COLS", "UNIT_UM")
                                                  for t in n.targets))]
    ns = {}
    exec(compile(ast.Module(body=keep, type_ignores=[]), "compare_runs.py", "exec"), ns)
    return ns


def load_spec(d, ns, col):
    pp = pd.read_csv(os.path.join(d, "cv_drug_ranks.tsv"), sep="\t")
    pp = ns["stratified_spec"](d, pp)
    pp = pp[pp.active.astype(str).str.upper() == "TRUE"]
    return pp.set_index(["fold", "pert"])[["drug", col]]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--qvals", required=True, help="expression limma qvals.csv (genes x perturbations, Ensembl ids)")
    ap.add_argument("--gem-genes", required=True, help="Human-GEM model/genes.tsv (column 'genes': Ensembl ids)")
    ap.add_argument("--run", action="append", required=True, metavar="NAME=DIR", help="real run; control: DIR_shuff")
    ap.add_argument("--q", type=float, default=0.05, help="DE threshold (with 0.05 the counts track the n_de_genes of the rank tables, r = 1.00)")
    ap.add_argument("--model", default="prob", choices=["prob", "logreg"])
    ap.add_argument("--stratum", default="timedose", choices=["time", "timedose"])
    ap.add_argument("--n-bins", type=int, default=4)
    ap.add_argument("--out-dir")
    args = ap.parse_args()
    runs = dict(r.split("=", 1) for r in args.run)
    out_dir = args.out_dir or os.path.join(os.path.dirname(list(runs.values())[0].rstrip("/")), "metabolic_coverage")
    os.makedirs(out_dir, exist_ok=True)
    col = f"spec_{args.model}_{args.stratum}_heldout"
    ns = compare_runs_functions()

    gem = set(pd.read_csv(args.gem_genes, sep="\t").genes)
    perts = set()
    spec = {}
    for name, d in runs.items():
        real, ctrl = load_spec(d, ns, col), load_spec(d.rstrip("/") + "_shuff", ns, col)
        j = real.join(ctrl[[col]], rsuffix="_ctrl", how="inner").dropna()
        spec[name] = j.rename(columns={col: "real", col + "_ctrl": "ctrl"})
        perts |= set(j.index.get_level_values("pert"))
    header = pd.read_csv(args.qvals, nrows=0).columns
    use = [header[0]] + [c for c in header[1:] if c in perts]
    q = pd.read_csv(args.qvals, index_col=0, usecols=use, engine="c")
    missing = perts - set(q.columns)
    if missing:
        print(f"{len(missing)} scored perturbations have no qvals column, left out (e.g. {sorted(missing)[:3]})")
    de = q < args.q
    in_gem = q.index.isin(gem)
    cov = pd.DataFrame({"n_de": de.sum(), "n_de_gem": de[in_gem].sum()})
    cov["frac_gem"] = cov.n_de_gem / cov.n_de.where(cov.n_de > 0)
    print(f"Human-GEM genes: {len(gem)}, in the qvals table: {in_gem.sum()} of {len(q)} "
          f"({in_gem.mean():.1%}, the share expected by chance among DE genes)")
    cov.rename_axis("pert").to_csv(os.path.join(out_dir, "coverage.tsv"), sep="\t")

    rows, rho = [], []
    for name, j in spec.items():
        j = j.join(cov, on="pert", how="inner")
        j["gain"] = j.ctrl - j.real  # positive: real better
        for var in ("n_de_gem", "frac_gem", "n_de"):
            r, p = spearmanr(j[var], j.gain, nan_policy="omit")
            rho.append(dict(run=name, variable=var, spearman_gain=r, p=p, n_perts=len(j)))
            bins = pd.qcut(j[var], args.n_bins, duplicates="drop")
            for b, g in j.groupby(bins, observed=True):
                dr = g.reset_index().groupby(["fold", g.reset_index().drug.str.lower()])[["real", "ctrl"]].mean()
                diff = dr.real - dr.ctrl
                rows.append(dict(run=name, binning=var, bin=str(b), n_perts=len(g), n_drugs=len(dr),
                                 real_pert=g.real.median(), ctrl_pert=g.ctrl.median(),
                                 real_drug=dr.real.median(), ctrl_drug=dr.ctrl.median(),
                                 gap_points=-100 * np.median(diff), frac_drugs_better=(diff < 0).mean(),
                                 wilcoxon_p=wilcoxon(diff).pvalue if (diff != 0).any() else np.nan))
    res = pd.DataFrame(rows)
    res.to_csv(os.path.join(out_dir, "spec_by_coverage.tsv"), sep="\t", index=False)
    pd.DataFrame(rho).to_csv(os.path.join(out_dir, "spearman.tsv"), sep="\t", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(res.round(4).to_string(index=False))
        print(pd.DataFrame(rho).round(4).to_string(index=False))
    print(f"written to {out_dir}")


if __name__ == "__main__":
    main()
