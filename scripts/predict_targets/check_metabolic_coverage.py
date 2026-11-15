#!/usr/bin/env python
"""
Does a signature's specificity depend on how much metabolic change a perturbation has?

Perturbations are called active from transcriptome-wide DE counts (--min-de-genes), but MetChange only sees the genes
of Human-GEM, and propagates them through the GEM into metabolite scores. This script measures metabolic coverage per
perturbation in up to two ways and splits the leak-free spec_<model>_<stratum>_heldout of each run and its shuffled
control by it (drugs as units: mean over the drug's perturbations in the bin; paired real - control gap; Wilcoxon p).

Coverage from the expression labels (--qvals, --gem-genes):
  n_de        DE genes (q < --q), transcriptome-wide
  n_de_gem    DE genes among the Human-GEM genes
  frac_gem    n_de_gem / n_de
Coverage from MetChange's own output (--mc-qvals; --mc-log2fc for the top-n columns; --mc-nodes to keep only the
features the model can use, i.e. those mapping to a metabolite node, as feature_node_map does):
  n_sig_mc    MetChange features with q < --q
  n_sig_topn  of the --top-n features with the largest |log2FC| (what the model reads with --ignore-qvalue
              --top-n), how many have q < --q. The model always reads --top-n features, significant or not; a
              perturbation with few significant ones gets a signature that is mostly noise.

  python check_metabolic_coverage.py --qvals <expr labels>/qvals.csv --gem-genes Human-GEM/model/genes.tsv \
      --mc-qvals <mc labels>/qvals.csv --mc-log2fc <mc labels>/log2FC.csv --mc-nodes <combined>/nodes.tsv \
      --run mc=<real out_dir> --run expr=<real out_dir> [--out-dir DIR]

Use the labels the runs used (on the cluster: the scratch copies). Each --run is a real run; its control is
<out_dir>_shuff. The spec percentiles are computed with compare_runs.py's own functions (same strata, mid-rank
percentiles), so they match compare_real_vs_shuffled/. Writes coverage.tsv (per perturbation), spec_by_coverage.tsv
(per run, coverage variable and bin) and spearman.tsv (per-perturbation gain vs each coverage variable).
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


def read_columns(path, perts):
    """features x perturbations CSV, only the columns of perts (first column: feature ids)."""
    header = pd.read_csv(path, nrows=0).columns
    use = [header[0]] + [c for c in header[1:] if c in perts]
    m = pd.read_csv(path, index_col=0, usecols=use, engine="c")
    m.index = m.index.astype(str)
    missing = perts - set(m.columns)
    if missing:
        print(f"{os.path.basename(os.path.dirname(path))}/{os.path.basename(path)}: {len(missing)} of {len(perts)} "
              f"scored perturbations have no column, left out (e.g. {sorted(missing)[:3]})")
    return m


def expression_coverage(args, perts):
    q = read_columns(args.qvals, perts)
    gem = set(pd.read_csv(args.gem_genes, sep="\t").genes)
    de = q < args.q
    in_gem = q.index.isin(gem)
    print(f"Human-GEM genes: {len(gem)}, in the expression labels: {in_gem.sum()} of {len(q)} ({in_gem.mean():.1%}, "
          f"the share of DE genes expected by chance)")
    cov = pd.DataFrame({"n_de": de.sum(), "n_de_gem": de[in_gem].sum()})
    cov["frac_gem"] = cov.n_de_gem / cov.n_de.where(cov.n_de > 0)
    return cov


def metchange_coverage(args, perts):
    q = read_columns(args.mc_qvals, perts)
    if args.mc_nodes:  # features the model can use: a metabolite node id, or a compartment-free id of one
        nodes = pd.read_csv(args.mc_nodes, sep="\t", dtype=str, keep_default_na=False)
        met = nodes[nodes.node_type == "metabolite"]
        usable = set(met.node_id) | set(met.base_met[met.base_met != ""])
        keep = q.index.isin(usable)
        print(f"MetChange features: {len(q)}, {keep.sum()} map to metabolite nodes (kept)")
        q = q[keep]
    sig = q < args.q
    cov = pd.DataFrame({"n_sig_mc": sig.sum()})
    if args.mc_log2fc:
        lfc = read_columns(args.mc_log2fc, perts).reindex(index=q.index, columns=q.columns)
        top = lfc.abs().rank(ascending=False, method="first") <= args.top_n
        cov["n_sig_topn"] = (sig & top).sum()
    return cov


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--qvals", help="expression limma qvals.csv (genes x perturbations, Ensembl ids)")
    ap.add_argument("--gem-genes", help="Human-GEM model/genes.tsv (column 'genes': Ensembl ids); with --qvals")
    ap.add_argument("--mc-qvals", help="MetChange limma qvals.csv (features x perturbations)")
    ap.add_argument("--mc-log2fc", help="MetChange log2FC.csv, for n_sig_topn; with --mc-qvals")
    ap.add_argument("--mc-nodes", help="network nodes.tsv: keep only MetChange features mapping to metabolite nodes")
    ap.add_argument("--top-n", type=int, default=100, help="as the runs' --top-n")
    ap.add_argument("--run", action="append", required=True, metavar="NAME=DIR", help="real run; control: DIR_shuff")
    ap.add_argument("--q", type=float, default=0.05,
                    help="significance threshold (expression: with 0.05 the counts track the n_de_genes of the rank "
                         "tables, r = 1.00)")
    ap.add_argument("--model", default="prob", choices=["prob", "logreg"])
    ap.add_argument("--stratum", default="timedose", choices=["time", "timedose"])
    ap.add_argument("--n-bins", type=int, default=4)
    ap.add_argument("--out-dir")
    args = ap.parse_args()
    if not args.qvals and not args.mc_qvals:
        ap.error("give --qvals (with --gem-genes) and/or --mc-qvals")
    if args.qvals and not args.gem_genes:
        ap.error("--qvals needs --gem-genes")
    if (args.mc_log2fc or args.mc_nodes) and not args.mc_qvals:
        ap.error("--mc-log2fc and --mc-nodes need --mc-qvals")
    runs = dict(r.split("=", 1) for r in args.run)
    out_dir = args.out_dir or os.path.join(os.path.dirname(list(runs.values())[0].rstrip("/")), "metabolic_coverage")
    os.makedirs(out_dir, exist_ok=True)
    col = f"spec_{args.model}_{args.stratum}_heldout"
    ns = compare_runs_functions()

    spec, perts = {}, set()
    for name, d in runs.items():
        real, ctrl = load_spec(d, ns, col), load_spec(d.rstrip("/") + "_shuff", ns, col)
        j = real.join(ctrl[[col]], rsuffix="_ctrl", how="inner").dropna()
        spec[name] = j.rename(columns={col: "real", col + "_ctrl": "ctrl"})
        perts |= set(j.index.get_level_values("pert"))

    parts = []
    if args.qvals:
        parts.append(expression_coverage(args, perts))
    if args.mc_qvals:
        parts.append(metchange_coverage(args, perts))
    cov = pd.concat(parts, axis=1, join="inner")
    variables = list(cov.columns)
    cov.rename_axis("pert").to_csv(os.path.join(out_dir, "coverage.tsv"), sep="\t")
    with pd.option_context("display.width", 200):
        print(cov.describe(percentiles=[.01, .05, .25, .5, .75, .95]).round(3).to_string())
        print("Spearman between coverage variables:")
        print(cov.corr("spearman").round(3).to_string())

    rows, rho = [], []
    for name, j in spec.items():
        j = j.join(cov, on="pert", how="inner")
        j["gain"] = j.ctrl - j.real  # positive: real better
        for var in variables:
            r, p = spearmanr(j[var], j.gain, nan_policy="omit")
            rho.append(dict(run=name, variable=var, spearman_gain=r, p=p, n_perts=len(j)))
            bins = pd.qcut(j[var], args.n_bins, duplicates="drop")
            for b, g in j.groupby(bins, observed=True):
                g = g.reset_index()
                dr = g.groupby(["fold", g.drug.str.lower()])[["real", "ctrl"]].mean()
                diff = dr.real - dr.ctrl
                rows.append(dict(run=name, variable=var, bin=str(b), n_perts=len(g), n_drugs=len(dr),
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
