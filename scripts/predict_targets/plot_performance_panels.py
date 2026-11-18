#!/usr/bin/env python3
"""Performance figure for the current recipe: MetChange, expr and TF against their own shuffled controls.

Four panels (reports/2026-11-18_recipe_state.md):
  A  ECDF of the percentile per drug, real vs own shuffled control, one panel per signature.
  B  Forest of the per-drug gap (control - real, points) with bootstrap CI, one row per run.
  C  Top-k recovery (true drug / target in the top 1 / 5 / 10 / 25 %), real bars with the control marked.
All panels are the same runs (whatever --mc/--expr/--tf name), except the forest, which also carries the
comparison rows given with --extra.

Gaps are per drug: each drug's percentile is averaged over its active perturbations in the run and in the run's
own `_shuff` control, then gap = 100 * (control - real). Lower percentile is better, so a positive gap means the
real signature ranks the drug's targets higher than the shuffled one. `spec_prob_timedose_heldout` by default:
the drug's own query against other queries of the same time x relative dose stratum, so the design cannot carry
drug identity.

Usage (from scripts/predict_targets, laptop or cluster):
    python plot_performance_panels.py --runs-dir ../../results/predict_targets/A549 \
        --coverage-dir ../../results/predict_targets/A549/metabolic_coverage_fixedlabels \
        --out-dir ../../results/predict_targets/A549/performance_figure
Needs each run's `cv_drug_ranks.tsv` and `fold*/test_scores_heldout.npz` (for the stratified spec metric), and
its `_shuff` counterpart. Runs missing a fold are used as they are, with the fold count reported in the caption
file; pass --require-folds 5 to skip incomplete runs instead.
"""
import argparse
import glob
import json
import os
import sys

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# Palette: slots 1-3 of the documented categorical theme (blue / orange / aqua), which the palette reference
# certifies for all-pairs use in both modes; chrome and ink from the same reference.
SERIES = {"mc": "#2a78d6", "expr": "#eb6834", "tf": "#1baf7a"}
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
AXIS = "#c3c2b7"
CONTROL = "#898781"  # shuffled control: muted, never a series hue

LABEL = {"mc": "MetChange", "expr": "expression", "tf": "TF"}


def metric_words(metric):
    """What the percentile is of, for the axis labels: spec_* ranks the drug's targets, rank_* the drug."""
    if metric.startswith("spec_"):
        return {"x": "percentile of the drug's own query (lower is better)",
                "head": "the drug's targets under its own query vs other queries, real vs own shuffled control",
                "c_y": "queries scoring their drug there (%)", "c_x": "specificity cut-off"}
    return {"x": "percentile of the true drug (lower is better)",
            "head": "rank of the true drug among candidates, real vs own shuffled control",
            "c_y": "drugs ranked there (%)", "c_x": "candidate ranking cut-off"}

#: dark steps of the same eight hues and the dark chrome / ink, selected for the dark surface (not a flip of
#: the light values): for a figure that will sit on a dark slide.
DARK = {"SERIES": {"mc": "#3987e5", "expr": "#d95926", "tf": "#199e70"}, "SURFACE": "#1a1a19",
        "INK": "#f0efec", "INK2": "#c3c2b7", "MUTED": "#898781", "GRID": "#2c2c2a", "AXIS": "#383835",
        "CONTROL": "#898781"}


def apply_theme(name):
    """Switch the module's colour roles to the light or dark set, then restate the rcParams from them."""
    if name == "dark":
        globals().update(DARK)
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "font.family": "sans-serif", "font.size": 9,
        "axes.edgecolor": AXIS, "axes.linewidth": 0.8, "axes.labelcolor": INK2, "axes.titlecolor": INK,
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "grid.color": GRID, "grid.linewidth": 0.6, "legend.frameon": False, "legend.labelcolor": INK2,
    })


apply_theme("light")


def load_ranks(run_dir):
    """cv_drug_ranks.tsv, or the per-fold files concatenated when the CV summary is missing (incomplete run)."""
    f = os.path.join(run_dir, "cv_drug_ranks.tsv")
    if os.path.exists(f):
        return pd.read_csv(f, sep="\t")
    parts = []
    for p in sorted(glob.glob(os.path.join(run_dir, "fold*", "test_drug_ranks.tsv"))):
        t = pd.read_csv(p, sep="\t")
        t["fold"] = int(os.path.basename(os.path.dirname(p)).replace("fold", ""))
        parts.append(t)
    if not parts:
        raise FileNotFoundError(f"{run_dir}: no cv_drug_ranks.tsv and no fold*/test_drug_ranks.tsv")
    return pd.concat(parts, ignore_index=True)


def per_drug(run_dir, stratified_spec, col):
    """Percentile of the true target per drug (mean over the drug's active perturbations).

    The stratified `spec_*` columns are already percentiles; the `rank_*` columns are raw ranks and are divided
    by the matching candidate count, exactly as `compare_runs.py` does (`n_spec*` / `n_candidates*`)."""
    pp = stratified_spec(run_dir, load_ranks(run_dir))
    if col not in pp.columns:
        raise SystemExit(f"{run_dir}: {col} not available (test_scores_heldout.npz missing in some fold?)")
    pp = pp[pp.active.astype(str).str.upper() == "TRUE"].copy()
    if col.startswith("spec_") and col.endswith("_heldout") and f"n_spec_heldout" in pp.columns and \
            pp[col].max(skipna=True) is not None and pp[col].max(skipna=True) > 1:
        pp["_pct"] = pp[col] / pp["n_spec_heldout"]          # raw spec rank (non-stratified column)
    elif col.startswith("rank_"):
        suffix = col[col.index("_heldout"):] if col.endswith("_heldout") else ""
        denom = "n_candidates" + suffix
        if denom not in pp.columns:
            raise SystemExit(f"{run_dir}: {denom} missing, cannot turn {col} into a percentile")
        pp["_pct"] = pp[col] / pp[denom]
    else:
        pp["_pct"] = pp[col]                                  # stratified spec_*: already a percentile
    pp = pp.dropna(subset=["_pct"])
    pp = pp.assign(drug=pp.drug.astype(str).str.lower())
    return pp.groupby("drug")["_pct"].mean(), sorted(pp.fold.unique())


def gaps(runs_dir, run, stratified_spec, col):
    """Per-drug gap in points against the run's own shuffled control, plus the folds present."""
    real, folds = per_drug(os.path.join(runs_dir, run), stratified_spec, col)
    ctrl, _ = per_drug(os.path.join(runs_dir, run + "_shuff"), stratified_spec, col)
    j = pd.concat([real.rename("real"), ctrl.rename("ctrl")], axis=1, join="inner").dropna()
    return j, 100.0 * (j.ctrl - j.real), folds


def boot_ci(x, n=5000, seed=0):
    rng = np.random.default_rng(seed)
    x = np.asarray(x, float)
    med = np.median(rng.choice(x, (n, len(x)), replace=True), axis=1)
    return np.median(x), np.percentile(med, 2.5), np.percentile(med, 97.5)


# ---------------------------------------------------------------- panels
def panel_a(axes, data, words):
    """ECDF of the true target's percentile: real vs own shuffled control, one axis per signature."""
    for ax, sig in zip(axes, ["mc", "expr", "tf"]):
        j = data[sig]["j"]
        ax.plot([0, 1], [0, 1], color=AXIS, lw=0.8, ls=(0, (4, 3)), zorder=1)
        for name, vals, color, lw in [("shuffled control", j.ctrl, CONTROL, 1.6),
                                      (f"{LABEL[sig]} signature", j.real, SERIES[sig], 2.0)]:
            v = np.sort(np.asarray(vals, float))
            ax.plot(v, np.arange(1, len(v) + 1) / len(v), color=color, lw=lw,
                    label=name, solid_capstyle="round", zorder=3)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel(words["x"], labelpad=4)
        ax.set_title(f"{LABEL[sig]}  (n = {len(j)} drugs)", fontsize=9.5, pad=6)
        ax.grid(axis="both", zorder=0)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.legend(loc="lower right", fontsize=7.5, handlelength=1.6)
    axes[0].set_ylabel("cumulative fraction of drugs")
    for ax in axes[1:]:
        ax.set_yticklabels([])


def panel_b(ax, rows):
    """Forest of the per-drug gap with bootstrap CI, grouped by signature."""
    y = np.arange(len(rows))[::-1]
    for yi, r in zip(y, rows):
        c = SERIES[r["sig"]]
        ax.plot([r["lo"], r["hi"]], [yi, yi], color=c, lw=2, solid_capstyle="round", zorder=3)
        ax.plot([r["med"]], [yi], "o", ms=6, color=c, mec=SURFACE, mew=1.2, zorder=4)
        ax.annotate(f"{r['med']:.1f}", (r["hi"], yi), xytext=(5, 0), textcoords="offset points",
                    va="center", fontsize=7.5, color=INK2)
    ax.axvline(0, color=AXIS, lw=0.8, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([r["label"] for r in rows], fontsize=8)
    ax.set_xlabel("per-drug gap over own shuffled control (points)")
    ax.grid(axis="x", zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="y", length=0)
    # no legend: each row is named on the y axis, so colour is redundant encoding here
    ax.set_xlim(left=min(0, min(r["lo"] for r in rows) - 0.5))
    ax.set_ylim(-0.9, len(rows) - 0.1)   # keep the rows together instead of spread over the whole box


def panel_c(ax, data, words, ks=(1, 5, 10, 25)):
    """Top-k recovery: fraction of drugs with the true target in the top k % of candidates."""
    sigs = ["mc", "expr", "tf"]
    w = 0.26
    x = np.arange(len(ks))
    for i, sig in enumerate(sigs):
        j = data[sig]["j"]
        off = (i - 1) * w
        real = [100 * np.mean(j.real <= k / 100) for k in ks]
        ctrl = [100 * np.mean(j.ctrl <= k / 100) for k in ks]
        ax.bar(x + off, real, w * 0.9, color=SERIES[sig], label=LABEL[sig], zorder=3,
               edgecolor=SURFACE, linewidth=1.2)
        for xi, c in zip(x + off, ctrl):
            ax.plot([xi - w * 0.33, xi + w * 0.33], [c, c], color=CONTROL, lw=2.2,
                    solid_capstyle="butt", zorder=5)
    ax.set_xticks(x)
    ax.set_xticklabels([f"top {k}%" for k in ks])
    ax.set_ylabel(words["c_y"])
    ax.set_xlabel(words["c_x"])
    ax.grid(axis="y", zorder=0)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    h = [Patch(facecolor=SERIES[s], label=LABEL[s]) for s in sigs]
    h.append(Line2D([0], [0], color=CONTROL, lw=2, label="own shuffled control"))
    ax.legend(handles=h, loc="upper left", fontsize=7.5)


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--runs-dir", required=True, help="cell directory holding the runs and their _shuff controls")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--metric", default="spec_prob_timedose_heldout")
    ap.add_argument("--mc", default="metchange_rank_noGlobals_lr1e-3_min15_qsigtop")
    ap.add_argument("--expr", default="expr_noGlobals_lr1e-3_min15_qsigtop")
    ap.add_argument("--tf", default="tf_active200_noGlobals_lr1e-3_min15_qsigtop")
    ap.add_argument("--extra", action="append", default=[], metavar="SIG=RUN:LABEL",
                    help="extra row for the forest, e.g. mc=metchange_rank_noGlobals_lr1e-3_min15:MetChange _min15 s0")
    ap.add_argument("--require-folds", type=int, default=0, help="skip runs with fewer folds than this")
    ap.add_argument("--background", choices=("transparent", "surface"), default="transparent",
                    help="transparent (default; for slides) or the theme's chart surface")
    ap.add_argument("--theme", choices=("light", "dark"), default="light",
                    help="ink and hues for a light slide (default) or a dark one; transparent either way")
    args = ap.parse_args()

    apply_theme(args.theme)
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import check_metabolic_coverage as cmc
    stratified_spec = cmc.compare_runs_functions()["stratified_spec"]

    os.makedirs(args.out_dir, exist_ok=True)
    main_runs = {"mc": args.mc, "expr": args.expr, "tf": args.tf}

    data, notes = {}, {}
    for sig, run in main_runs.items():
        j, g, folds = gaps(args.runs_dir, run, stratified_spec, args.metric)
        if args.require_folds and len(folds) < args.require_folds:
            raise SystemExit(f"{run}: only folds {folds}, --require-folds {args.require_folds}")
        tag = "_" + run.split("_min15")[-1].lstrip("_") if run.split("_min15")[-1].strip("_") else "_min15"
        data[sig] = {"run": run, "j": j, "gap": g, "folds": folds, "tag": tag}
        notes[run] = {"signature": sig, "n_drugs": int(len(j)), "folds": [int(f) for f in folds],
                      "median_pct_real": float(j.real.median()), "median_pct_control": float(j.ctrl.median()),
                      "median_gap_points": float(g.median()),
                      "frac_drugs_better": float((g > 0).mean())}
        print(f"{run}: folds {folds}, {len(j)} drugs, gap {g.median():.1f} points")

    rows = []
    for sig in ("mc", "expr", "tf"):
        d = data[sig]
        med, lo, hi = boot_ci(d["gap"])
        tag = "" if len(d["folds"]) == 5 else f", folds {''.join(str(f) for f in d['folds'])}"
        rows.append({"sig": sig, "label": f"{LABEL[sig]}  (n = {len(d['gap'])}{tag})",
                     "med": med, "lo": lo, "hi": hi})
        notes[d["run"]].update(ci_low=lo, ci_high=hi)
    for spec in args.extra:
        sig, rest = spec.split("=", 1)
        run, label = rest.split(":", 1)
        j, g, folds = gaps(args.runs_dir, run, stratified_spec, args.metric)
        med, lo, hi = boot_ci(g)
        rows.append({"sig": sig, "label": f"{label} (n={len(g)})", "med": med, "lo": lo, "hi": hi})
        notes[run] = {"signature": sig, "n_drugs": int(len(j)), "folds": [int(f) for f in folds],
                      "median_gap_points": float(med), "ci_low": lo, "ci_high": hi,
                      "frac_drugs_better": float((g > 0).mean())}
        print(f"{run}: folds {folds}, {len(g)} drugs, gap {med:.1f} points")
    rows.sort(key=lambda r: (["mc", "expr", "tf"].index(r["sig"]), -r["med"]))

    fig = plt.figure(figsize=(11, 5.95))
    gs = fig.add_gridspec(2, 6, height_ratios=[1, 0.78], hspace=0.34, wspace=0.8,
                          left=0.075, right=0.985, top=0.925, bottom=0.10)
    axes_a = [fig.add_subplot(gs[0, 2 * i:2 * i + 2]) for i in range(3)]
    words = metric_words(args.metric)
    panel_a(axes_a, data, words)
    ax_b = fig.add_subplot(gs[1, 0:3])
    panel_b(ax_b, rows)
    ax_c = fig.add_subplot(gs[1, 3:6])
    panel_c(ax_c, data, words)

    axes_a[1].annotate(words["head"],
                       xy=(0.5, 1), xycoords="axes fraction", xytext=(0, 30), textcoords="offset points",
                       ha="center", fontsize=10.5, color=INK, annotation_clip=False)

    stem = "performance_panels_" + args.metric.replace("_heldout", "") + ("_dark" if args.theme == "dark" else "")
    for ext in ("png", "pdf"):
        p = os.path.join(args.out_dir, f"{stem}.{ext}")
        fig.savefig(p, dpi=300, bbox_inches="tight", transparent=args.background == "transparent")
        print("saved", p)
    with open(os.path.join(args.out_dir, f"{stem}.json"), "w") as fh:
        json.dump({"metric": args.metric, "runs": notes}, fh, indent=2)

    tbl = pd.DataFrame([{"run": k, **v} for k, v in notes.items()])
    tbl.to_csv(os.path.join(args.out_dir, f"{stem}.tsv"), sep="\t", index=False)
    print("saved", os.path.join(args.out_dir, f"{stem}.tsv"))


if __name__ == "__main__":
    main()
