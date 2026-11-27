#!/usr/bin/env python3
"""Effective dimensionality of a label matrix, and set-specific signal within it.

Report 11, section 8. Answers three questions:

  1. How many independent dimensions does a label matrix actually have?
  2. Is that number unusual, or is it what the input expression looks like too?
  3. Do the curated metabolite sets carry signal beyond the dominant global
     component, relative to size-matched random sets?

Written for interactive use: every step is a function with working defaults, so

    from feature_dimensionality import *
    X, ids, cols = load_matrix(label_path('HEPG2', 'metchange'))
    spectrum(X)

works line by line in Spyder or a VS Code `# %%` cell. `main()` only wires the
same functions to argparse.

WHY THE OBVIOUS PCA IS WRONG
----------------------------
Centring the standardised matrix across FEATURES (`Z - Z.mean(0)`) subtracts,
for each perturbation, its mean over all metabolites. When the features share a
dominant component that mean IS the shared component, so centring removes the
very thing being measured and PC1 drops from 95% to ~20%. This module always
works from the feature x feature correlation matrix, which has no such step.
That bug produced a wrong figure (9.55 effective dimensions instead of 1.10)
in the first pass at this analysis; the guard is `_check_no_double_centring`.
"""

import argparse
import csv
import os

import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SETS_DEFAULT = os.path.join(
    REPO, "results/mc_phenotype_assoc/metabolite_sets/metabolite_sets.tsv"
)
# The sets whose results reports 09 and 10 rest on.
SETS_OF_INTEREST = [
    ("dntp_endpoints", "all"),
    ("ntp_endpoints", "all"),
    ("nmp_endpoints", "all"),
    ("tca_endpoints", "all_nohub"),
    ("amino_acid_endpoints", "all"),
    ("oxphos", "mito"),
    ("glycolysis", "cyto_nohub"),
    ("etc_endpoints", "all_nohub"),
]


# On the cluster the label matrices live on scratch, not in the repo; on a laptop they
# are in the repo. `mc_phenotype_assoc_common.sh` resolves this for the shell launchers
# as `labels_base`, and reproduce_report11.sh exports it as LABELS_BASE for this module.
LABELS_BASE = os.environ.get("LABELS_BASE") or REPO


def label_path(cell="HEPG2", arm="metchange", stat="log2FC", base=None):
    """Path to one label matrix. arm is 'metchange' or 'expr'.

    Honours $LABELS_BASE so the same code runs against the repo (laptop) or scratch
    (cluster) without editing. Pass `base` to override explicitly.
    """
    sub = "metchange/rank" if arm == "metchange" else "expr"
    return os.path.join(base or LABELS_BASE, "results/labels", cell, sub, "%s.csv" % stat)


def load_matrix(path, max_rows=None, min_finite=0.9):
    """Read a label matrix. Returns (X, feature_ids, perturbation_names).

    Rows with more than (1 - min_finite) missing are dropped, then any column
    with a remaining NaN, so X is complete and no imputation is implied.
    """
    X, ids = [], []
    with open(path) as fh:
        header = next(fh).rstrip("\n").split(",")
        for line in fh:
            q = line.rstrip("\n").split(",")
            v = np.array(
                [np.nan if x.strip() in ("", "NA", "nan") else float(x) for x in q[1:]],
                dtype=float,
            )
            if np.isfinite(v).sum() > min_finite * v.size:
                ids.append(q[0])
                X.append(v)
            if max_rows and len(X) >= max_rows:
                break
    X = np.vstack(X)
    keep = np.isfinite(X).all(0)
    cols = [h for h, k in zip(header[1:], keep) if k]
    return X[:, keep], ids, cols


def standardise(X):
    """z-score each FEATURE across perturbations. No centring across features."""
    return (X - X.mean(1, keepdims=True)) / (X.std(1, keepdims=True) + 1e-12)


def _check_no_double_centring(Z):
    """The shared component must still be present; see the module docstring."""
    resid = np.abs(Z.mean(0)).mean()
    if resid < 1e-8:
        raise ValueError(
            "matrix appears centred across features: the shared component has "
            "already been removed, so PC1 will be understated"
        )


def spectrum(X, n_features=1500, seed=0):
    """Eigen-spectrum of the feature x feature correlation matrix.

    Returns a dict with pc1, pc1_3, n_eig_90, participation_ratio, median_r,
    frac_r_above_95. participation_ratio = (sum L)^2 / sum(L^2) is a continuous
    measure of how many dimensions carry the variance; it uses the whole
    spectrum and applies no variance cutoff.
    """
    Z = standardise(X)
    _check_no_double_centring(Z)
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(Z), min(n_features, len(Z)), replace=False)
    C = np.corrcoef(Z[idx])
    ev = np.clip(np.linalg.eigvalsh(C)[::-1], 0, None)
    var = ev / ev.sum()
    off = C[np.triu_indices(len(idx), 1)]
    return {
        "n_features": X.shape[0],
        "n_perturbations": X.shape[1],
        "n_sampled": len(idx),
        "pc1": 100 * var[0],
        "pc1_3": 100 * var[:3].sum(),
        "n_eig_90": int(np.searchsorted(np.cumsum(var), 0.90)) + 1,
        "participation_ratio": ev.sum() ** 2 / (ev ** 2).sum(),
        "median_r": float(np.median(off)),
        "frac_r_above_95": float(np.mean(off > 0.95)),
    }


def global_component(X):
    """The dominant shared direction: mean standardised score over all features."""
    Z = standardise(X)
    g = Z.mean(0)
    return (g - g.mean()) / g.std()


def variance_explained_by_global(X, stride=100):
    """Mean fraction of a single feature's variance explained by g alone."""
    Z = standardise(X)
    g = global_component(X)
    return float(np.mean([np.corrcoef(z, g)[0, 1] ** 2 for z in Z[::stride]]))


def read_sets(path=SETS_DEFAULT):
    """{(set, variant): [metabolite, ...]} from metabolite_sets.tsv."""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            out[(r["set"], r["variant"])] = r["members"].split(",")
    return out


def set_specificity(X, ids, members, n_random=2000, seed=0):
    """How much of a set's score is NOT the global component, versus random sets.

    The set score is the mean standardised score over members, itself
    standardised. Regressing out g leaves a residual whose SD measures
    set-specific movement. A size-matched random set averages over unrelated
    metabolites and so collapses toward g, giving a smaller residual; a
    functionally coherent set has correlated deviations that survive averaging.

    A HIGH percentile therefore means the set carries real set-specific signal.
    It does NOT require the set to be orthogonal to g, which would be the wrong
    criterion: if the biology is partly aligned with the dominant direction,
    demanding orthogonality discards the biology.
    """
    pos = {m: i for i, m in enumerate(ids)}
    rows = [pos[m] for m in members if m in pos]
    if len(rows) < 2:
        return None
    Z = standardise(X)
    g = global_component(X)

    def residual_sd(block):
        s = block.mean(0)
        s = (s - s.mean()) / s.std()
        r = np.corrcoef(s, g)[0, 1]
        return r, (s - r * g).std()

    r_set, sd_set = residual_sd(Z[rows])
    rng = np.random.default_rng(seed)
    null = np.empty(n_random)
    for k in range(n_random):
        ix = rng.choice(len(Z), len(rows), replace=False)
        null[k] = residual_sd(Z[ix])[1]
    return {
        "n_members": len(rows),
        "n_requested": len(members),
        "r_with_global": float(r_set),
        "residual_sd": float(sd_set),
        "rand_percentile": 100 * float(np.mean(null < sd_set)),
    }


def compare_arms(cells=("HEPG2", "A549"), expr_max_rows=8000):
    """Question 1 and 2: is the dimensionality unusual for this data?"""
    print("effective dimensionality, from the feature x feature correlation matrix")
    print(
        "  %-30s %14s %7s %7s %9s %9s"
        % ("matrix", "shape", "PC1", "PR", "eig@90%", "median r")
    )
    out = {}
    for cell in cells:
        for arm, cap in (("metchange", None), ("expr", expr_max_rows)):
            path = label_path(cell, arm)
            if not os.path.exists(path):
                continue
            X, _, _ = load_matrix(path, max_rows=cap)
            s = spectrum(X)
            out[(cell, arm)] = s
            print(
                "  %-30s %6dx%-7d %6.1f%% %7.2f %9d %9.4f"
                % (
                    "%s  %s" % (cell, arm),
                    s["n_features"],
                    s["n_perturbations"],
                    s["pc1"],
                    s["participation_ratio"],
                    s["n_eig_90"],
                    s["median_r"],
                )
            )
    rng = np.random.default_rng(1)
    shape = (1500, 5000)
    s = spectrum(rng.standard_normal(shape))
    print(
        "  %-30s %6dx%-7d %6.1f%% %7.2f %9d %9.4f"
        % ("iid Gaussian (calibration)", shape[0], shape[1], s["pc1"],
           s["participation_ratio"], s["n_eig_90"], s["median_r"])
    )
    return out


def report_sets(cell="HEPG2", sets_path=SETS_DEFAULT, keys=None, n_random=2000):
    """Question 3: do the curated sets beat size-matched random sets?"""
    X, ids, _ = load_matrix(label_path(cell, "metchange"))
    sets = read_sets(sets_path)
    print()
    print("%s: variance of a single metabolite explained by the global component: %.1f%%"
          % (cell, 100 * variance_explained_by_global(X)))
    print()
    print("  %-26s %5s %10s %10s %13s"
          % ("set/variant", "n", "r(set,g)", "resid SD", "rand pctile"))
    out = {}
    for key in keys or SETS_OF_INTEREST:
        res = set_specificity(X, ids, sets.get(key, []), n_random=n_random)
        if res is None:
            print("  %-26s   (fewer than 2 members present)" % "/".join(key))
            continue
        out[key] = res
        print("  %-26s %5d %10.5f %10.4f %12.1f%%"
              % ("/".join(key), res["n_members"], res["r_with_global"],
                 res["residual_sd"], res["rand_percentile"]))
    print()
    print("  resid SD    = part of the set score that is NOT the global component")
    print("  rand pctile = its rank among size-matched random metabolite sets;")
    print("                HIGH means genuine set-specific signal (see docstring)")
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cells", nargs="+", default=["HEPG2", "A549"],
                    help="cell lines with local label matrices")
    ap.add_argument("--sets", default=SETS_DEFAULT)
    ap.add_argument("--n-random", type=int, default=2000,
                    help="size-matched random sets for the specificity null")
    ap.add_argument("--skip-sets", action="store_true")
    a = ap.parse_args()
    compare_arms(cells=tuple(a.cells))
    if not a.skip_sets:
        report_sets(cell=a.cells[0], sets_path=a.sets, n_random=a.n_random)


if __name__ == "__main__":
    main()
