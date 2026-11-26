#!/usr/bin/env python3
"""Build the two derived label matrices report 11 sections 9 and 10 compare against.

Both are written as a drop-in `--labels` directory for `mc_phenotype_assoc.py`, so
every downstream analysis (set scoring, covariate matching, random-set nulls) runs
unchanged on them.

  --mode producer_average
      The metabolite-centric EXPRESSION comparator (report 11, section 6). For each
      metabolite, the mean expression log2FC of the genes catalysing reactions that
      can PRODUCE it. Same gene-level information MetChange starts from, no
      stoichiometry, no LP, no network propagation. If MetChange cannot beat this,
      the integration is adding nothing that averaging the relevant enzymes does not.

      This is the right comparator here, and E-Flux / GIMME / iMAT are not, because
      those are reaction-centric: they score fluxes, and the quantity under test is a
      per-metabolite production capacity. Validating MetChange as an FBA method was
      Zielinski's work and is already published; the question here is whether the
      integration models metabolic change well enough to justify it as the FNN's
      metabolic signature.

  --mode global_residual
      MetChange with the dominant global component regressed out PER FEATURE
      (report 11, section 8). Each feature is standardised, the global component g
      (the mean standardised score across all features) is computed once, and each
      feature's own loading on g is removed:  resid[m] = Z[m] - corr(Z[m], g) * g.

      This is NOT what the `z` statistic already does. `z` subtracts a per-column
      UNIFORM mean, identical for every feature; the loadings differ by feature, so
      `z` removes only the uniform part. Measured: set-level z computed on the raw
      versus residualised matrix correlates 0.26 (oxphos/mito) to 0.85
      (glycolysis/cyto_nohub). The existing results are partially, not fully,
      protected against the global component.

Units differ between the two and from the source: `producer_average` is in
expression log2FC units, `global_residual` in per-feature SD units. Every statistic
in `mc_phenotype_assoc.py` is either column-standardised (`z`) or rank-based
(`auroc`), so neither is comparable in absolute level but both are comparable in
the statistics actually reported.

Interactive use, like `feature_dimensionality.py`:

    from build_derived_labels import *
    prod = producer_genes()                       # metabolite -> {ENSG}
    M, mets, cols = producer_average_matrix('HEPG2', prod)
"""

import argparse
import csv
import os
import shutil
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(HERE))
GEM_XML = os.path.join(REPO, "data/GEMs/Human-GEM/model/Human-GEM.xml")

sys.path.insert(0, HERE)
from feature_dimensionality import (  # noqa: E402
    label_path,
    load_matrix,
    standardise,
    global_component,
)

# Carried over from the GEM-derived sets: these appear in thousands of reactions,
# so their "producers" would be most of the model and the producer average would
# degenerate to the global expression mean. Reported separately rather than dropped.
CURRENCY = {
    "MAM01371", "MAM01285", "MAM01334", "MAM02040", "MAM02039", "MAM02751",
    "MAM02759", "MAM01948", "MAM01802", "MAM02553", "MAM02552", "MAM01597",
    "MAM03103", "MAM02630", "MAM02631", "MAM01253", "MAM02877", "MAM02519",
    "MAM01688", "MAM02046", "MAM02819",
}


def producer_genes(gem_xml=GEM_XML, include_currency=False, verbose=True):
    """metabolite id -> {Ensembl gene id} for reactions that can PRODUCE it.

    A reaction produces metabolite m if m is on the product side and the reaction
    can carry forward flux, or on the substrate side and it can carry reverse flux.
    Reactions with no GPR (exchanges, diffusion, sinks) contribute no genes, which
    is correct here: they carry no expression information.
    """
    from cobra.io import read_sbml_model

    gem = read_sbml_model(gem_xml)
    out = {}
    n_norule = 0
    for rxn in gem.reactions:
        genes = {g.id for g in rxn.genes}
        if not genes:
            n_norule += 1
            continue
        fwd = rxn.upper_bound > 0
        rev = rxn.lower_bound < 0
        for met, coef in rxn.metabolites.items():
            produced = (coef > 0 and fwd) or (coef < 0 and rev)
            if not produced:
                continue
            if not include_currency and met.id[:8] in CURRENCY:
                continue
            out.setdefault(met.id, set()).update(genes)
    if verbose:
        sizes = np.array([len(v) for v in out.values()])
        print(
            "producers: %d metabolites with >=1 gene; median %d genes, "
            "90th pct %d, max %d; %d reactions had no GPR"
            % (len(out), np.median(sizes), np.percentile(sizes, 90), sizes.max(), n_norule),
            file=sys.stderr,
        )
    return out


def arm_columns(cell, arm):
    """Perturbation names of one label matrix, in file order, without loading it."""
    with open(label_path(cell, arm)) as fh:
        return fh.readline().rstrip("\n").split(",")[1:]


def shared_columns(cell, arms=("metchange", "expr"), verbose=True):
    """Perturbations present in EVERY arm, in the first arm's order.

    THE ARMS ARE NOT ALWAYS CO-EXTENSIVE AND THE MISMATCH CAN BE SEVERE. Measured
    2026-11-26: HEPG2 has the identical 15,032 perturbations in both arms, in the
    same order, so a comparison there needs no restriction. A549 has 44,396 in the
    MetChange arm against 7,636 in the expression arm, sharing 7,633 - so comparing
    the arms unrestricted would score MetChange on nearly six times as many
    perturbations, and hence a different drug set, as the expression arms. Any
    difference between arms would then partly be a difference in which drugs were
    included. Restrict first; `reproduce_report11.sh` does this automatically.
    """
    cols = [arm_columns(cell, a) for a in arms]
    keep = set(cols[0]).intersection(*[set(c) for c in cols[1:]])
    out = [c for c in cols[0] if c in keep]
    if verbose:
        print("shared columns for %s: %s -> %d"
              % (cell, " / ".join("%s %d" % (a, len(c)) for a, c in zip(arms, cols)),
                 len(out)), file=sys.stderr)
    return out


def restrict_matrix(cell, arm, columns, verbose=True):
    """(M, ids, cols): one label matrix cut down to `columns`, in that order."""
    X, ids, cols = load_matrix(label_path(cell, arm))
    ci = {c: i for i, c in enumerate(cols)}
    missing = [c for c in columns if c not in ci]
    take = [ci[c] for c in columns if c in ci]
    if verbose and missing:
        print("  %d of %d requested columns absent after the NaN filter (dropped)"
              % (len(missing), len(columns)), file=sys.stderr)
    kept = [c for c in columns if c in ci]
    if verbose:
        print("restricted %s/%s to %d perturbations" % (cell, arm, len(kept)),
              file=sys.stderr)
    return X[:, take], ids, kept


def producer_average_matrix(cell, producers=None, expr_max_rows=None, verbose=True,
                            columns=None):
    """(M, metabolite_ids, perturbations): mean expression log2FC over producers."""
    if producers is None:
        producers = producer_genes()
    if columns is None:
        E, genes, cols = load_matrix(label_path(cell, "expr"), max_rows=expr_max_rows)
    else:
        E, genes, cols = restrict_matrix(cell, "expr", columns)
    gi = {g: i for i, g in enumerate(genes)}
    # Restrict to the metabolites the MetChange matrix carries, so the metabolite
    # sets index both matrices identically.
    _, mc_ids, _ = load_matrix(label_path(cell, "metchange"), max_rows=1)
    with open(label_path(cell, "metchange")) as fh:
        next(fh)
        mc_ids = [line.split(",", 1)[0] for line in fh]
    rows, kept, dropped = [], [], 0
    for m in mc_ids:
        idx = [gi[g] for g in producers.get(m, ()) if g in gi]
        if not idx:
            dropped += 1
            continue
        rows.append(E[idx].mean(0))
        kept.append(m)
    if verbose:
        print(
            "producer-average: %d of %d metabolites scored (%d had no measured "
            "producer gene), %d perturbations"
            % (len(kept), len(mc_ids), dropped, len(cols)),
            file=sys.stderr,
        )
    return np.vstack(rows), kept, cols


def global_residual_matrix(cell, verbose=True, columns=None):
    """(M, metabolite_ids, perturbations): MetChange with g removed per feature.

    With `columns`, the global component is computed on the RESTRICTED matrix, which
    is the right thing: g is a property of the perturbation set being analysed.
    """
    if columns is None:
        X, ids, cols = load_matrix(label_path(cell, "metchange"))
    else:
        X, ids, cols = restrict_matrix(cell, "metchange", columns)
    Z = standardise(X)
    g = global_component(X)
    load = np.array([np.corrcoef(z, g)[0, 1] for z in Z])
    R = Z - load[:, None] * g[None, :]
    if verbose:
        print(
            "global residual: %d features x %d perturbations; loadings on g "
            "median %.4f [min %.4f, max %.4f]; variance remaining %.1f%%"
            % (R.shape[0], R.shape[1], np.median(load), load.min(), load.max(),
               100 * (R.var(1) / Z.var(1)).mean()),
            file=sys.stderr,
        )
    return R, ids, cols


def write_labels(M, ids, cols, out_dir, source_dir, verbose=True):
    """Write log2FC.csv and copy the auxiliary files mc_phenotype_assoc.py expects."""
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, "log2FC.csv")
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["V1"] + list(cols))
        for i, m in enumerate(ids):
            w.writerow([m] + ["%.7g" % v for v in M[i]])
    # qvals and the n_sig diagnostics belong to the ORIGINAL limma fit. They are
    # carried over ONLY when this matrix has the source's exact perturbation set;
    # otherwise they describe different columns and must not be present at all.
    #
    # This is not hypothetical. Carrying them over blindly broke the A549 run on
    # 2026-11-26: the restricted matrices had 7,633 columns while the symlinked
    # qvals had 44,397 (metchange) and 7,637 (expr). The expression arm died with
    # "matrix q has a different header to value"; the MetChange arm did NOT die,
    # because its symlinked fit_info_nsig.csv satisfied the n_sig lookup first - so
    # it silently took n_sig from a different perturbation set. A loud failure in
    # one arm and a silent inconsistency in another, from the same cause.
    src_cols = None
    src_matrix = os.path.join(source_dir, "log2FC.csv")
    if os.path.exists(src_matrix):
        with open(src_matrix) as fh:
            src_cols = fh.readline().rstrip("\n").split(",")[1:]
    same_cols = src_cols is not None and list(cols) == src_cols
    carried = []
    if same_cols:
        for aux in ("qvals.csv", "fit_info_nsig.csv", "fit_info.csv"):
            src = os.path.join(source_dir, aux)
            dst = os.path.join(out_dir, aux)
            if os.path.exists(src) and not os.path.exists(dst):
                try:
                    os.symlink(os.path.relpath(src, out_dir), dst)
                except OSError:
                    shutil.copy2(src, dst)
                carried.append(aux)
    else:
        for aux in ("qvals.csv", "fit_info_nsig.csv", "fit_info.csv"):
            dst = os.path.join(out_dir, aux)
            if os.path.islink(dst) or os.path.exists(dst):
                os.remove(dst)
    with open(os.path.join(out_dir, "PROVENANCE.txt"), "w") as fh:
        fh.write(
            "Derived matrix, not a limma fit. Built by build_derived_labels.py.\n"
            "source: %s\n%d features x %d perturbations\n\n" % (source_dir, M.shape[0], M.shape[1]))
        if same_cols:
            fh.write(
                "Perturbation set is IDENTICAL to the source, so these are symlinked\n"
                "from it: %s\nThey describe the original limma fit, NOT this matrix:\n"
                "use --mode value only.\n" % (", ".join(carried) or "none present"))
        else:
            fh.write(
                "Perturbation set DIFFERS from the source (%d here, %d there), so\n"
                "qvals.csv and fit_info*.csv are deliberately ABSENT: they index\n"
                "different columns and would be wrong. Run with --mode value and\n"
                "--no-nsig.\n" % (M.shape[1], len(src_cols) if src_cols else -1))
    if verbose:
        print("wrote %s (%d x %d)" % (path, M.shape[0], M.shape[1]), file=sys.stderr)
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cell", default="HEPG2")
    ap.add_argument("--mode", required=True,
                    choices=["producer_average", "global_residual", "metchange_shared"])
    ap.add_argument("--restrict-to-shared", action="store_true",
                    help="cut every arm down to the perturbations present in BOTH the "
                         "metchange and expr matrices. Needed whenever the arms are not "
                         "co-extensive (A549: 44,396 vs 7,636) or the comparison is not "
                         "like for like; a no-op where they already match (HEPG2).")
    ap.add_argument("--out", help="output label dir (default: results/labels/<cell>/<mode>)")
    ap.add_argument("--include-currency", action="store_true",
                    help="producer_average: keep ATP/NAD/H2O etc., whose producer "
                         "sets are most of the model")
    a = ap.parse_args()
    out = a.out or os.path.join(REPO, "results/labels", a.cell, a.mode)
    source = os.path.dirname(label_path(a.cell, "metchange"))
    cols_want = shared_columns(a.cell) if a.restrict_to_shared else None
    if a.mode == "producer_average":
        M, ids, cols = producer_average_matrix(
            a.cell, producer_genes(include_currency=a.include_currency),
            columns=cols_want)
    elif a.mode == "global_residual":
        M, ids, cols = global_residual_matrix(a.cell, columns=cols_want)
    else:
        if cols_want is None:
            cols_want = shared_columns(a.cell)
        M, ids, cols = restrict_matrix(a.cell, "metchange", cols_want)
    write_labels(M, ids, cols, out, source)


if __name__ == "__main__":
    main()
