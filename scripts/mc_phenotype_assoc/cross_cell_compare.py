#!/usr/bin/env python3
"""Compare the phenotype association across cell lines, on a COMMON drug set.

Why the restriction is not optional. Every cell line has labels for a different
subset of the mapped compounds - 717 drugs for HEPG2 against 1,364 for PC3 - so
the raw per-cell numbers confound coverage with biology. Restricting to the drugs
present in EVERY cell line with one consistent group makes the cell line the only
thing that varies. It matters in practice: PC3's MMP-down association is
nominally significant on its own drugs (AUROC 0.454, p = 0.036) and gone on the
common set (0.478, p = 0.46).

Two designs are reported per cell:
  pooled    Mann-Whitney on the common set.
  matched   the same after 1:1 greedy matching of control to case drugs on
            log1p(n_sig) within that cell line, because the responsive share
            differs widely between cell lines (case/control 28%/7% in HEPG2,
            64%/37% in HA1E) and report 08 section 3 showed responsiveness is the
            confounder that matters. HA1E's nominal MMP-down signal (p = 0.025)
            does not survive it; HEPG2's does.

Reads only per_drug_scores.tsv, so it is cheap and needs no label matrices.
Results of the first run: report 09 (reports/2026-11-23_cross_cell_specificity.md).
"""
import argparse, csv, os, sys

import numpy as np

try:
    from scipy import stats as sps
except Exception:                                      # noqa: BLE001
    sps = None

DEFAULT_CELLS = "HEPG2 A549 A375 HA1E HCC515 HT29 PC3 VCAP"


def load(base, cell, run, set_name, variant, statistic):
    p = os.path.join(base, cell, run, "per_drug_scores.tsv")
    if not os.path.isfile(p):
        return {}
    out = {}
    for r in csv.DictReader(open(p), delimiter="\t"):
        if r["set"] != set_name or r["variant"] != variant:
            continue
        try:
            v = float(r[statistic])
            ns = int(r["max_n_sig"])
        except (ValueError, KeyError):
            continue
        out[r["drug"]] = (r["group"], v, ns)
    return out


def greedy_match(case, ctrl, caliper):
    """1:1 nearest-neighbour match on log1p(n_sig), hardest case first, no replacement."""
    pool = list(ctrl)
    med = float(np.median([c[1] for c in ctrl])) if ctrl else 0.0
    pairs = []
    for k in sorted(range(len(case)), key=lambda k: -abs(case[k][1] - med)):
        if not pool:
            break
        d = [abs(c[1] - case[k][1]) for c in pool]
        j = int(np.argmin(d))
        if d[j] <= caliper:
            pairs.append((case[k][0], pool[j][0]))
            pool.pop(j)
    return pairs


def mw(a, b):
    a = np.asarray([x for x in a if np.isfinite(x)], float)
    b = np.asarray([x for x in b if np.isfinite(x)], float)
    if len(a) < 5 or len(b) < 5 or sps is None:
        return np.nan, np.nan, len(a), len(b)
    u, p = sps.mannwhitneyu(a, b, alternative="two-sided")
    return u / (len(a) * len(b)), p, len(a), len(b)


def main():
    repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--assoc-base", default=os.path.join(repo, "results/mc_phenotype_assoc"),
                    help="directory holding <cell>/<run>/per_drug_scores.tsv")
    ap.add_argument("--cells", default=DEFAULT_CELLS)
    ap.add_argument("--run", default="attene_allmodes")
    ap.add_argument("--set", dest="set_name", default="tca")
    ap.add_argument("--variant", default="mito_nohub")
    ap.add_argument("--statistic", default="z")
    ap.add_argument("--case", default="MMP_down")
    ap.add_argument("--control", default="inactive")
    ap.add_argument("--caliper", type=float, default=0.25)
    ap.add_argument("--out")
    args = ap.parse_args()

    cells = [c for c in args.cells.replace(",", " ").split() if c]
    D = {c: load(args.assoc_base, c, args.run, args.set_name, args.variant, args.statistic)
         for c in cells}
    missing = [c for c in cells if not D[c]]
    D = {c: v for c, v in D.items() if v}
    if not D:
        raise SystemExit("no per_drug_scores.tsv with set=%s variant=%s under %s"
                         % (args.set_name, args.variant, args.assoc_base))
    if missing:
        print("no results for: %s" % ", ".join(missing), file=sys.stderr)

    groups = {}
    for c in D:
        for d, (g, _, _) in D[c].items():
            groups.setdefault(d, set()).add(g)
    common = [d for d in groups if all(d in D[c] for c in D) and len(groups[d]) == 1]
    cc = [d for d in common if next(iter(groups[d])) == args.case]
    ck = [d for d in common if next(iter(groups[d])) == args.control]

    print("%s %s / %s / %s vs %s" % (args.set_name, args.variant, args.statistic,
                                     args.case, args.control))
    print("drugs per cell: %s" % {c: len(D[c]) for c in D})
    print("common to all %d cells with one consistent group: %d drugs (%d case, %d control)"
          % (len(D), len(common), len(cc), len(ck)))
    print()
    hdr = "%-9s %22s %26s %18s" % ("cell", "pooled AUROC (p)", "matched AUROC (p, pairs)",
                                   "responsive c/ctrl")
    print(hdr)
    print("-" * len(hdr))
    rows = []
    for c in cells:
        if c not in D:
            continue
        a = [D[c][d][1] for d in cc]
        b = [D[c][d][1] for d in ck]
        pa, pp, na, nb = mw(a, b)
        ns = lambda ds: [(d, float(np.log1p(D[c][d][2]))) for d in ds
                         if np.isfinite(D[c][d][1])]
        pairs = greedy_match(ns(cc), ns(ck), args.caliper)
        if len(pairs) >= 20:
            ma, mp, _, _ = mw([D[c][x][1] for x, _ in pairs],
                              [D[c][y][1] for _, y in pairs])
            mtxt = "%.3f (p=%.2g, n=%d)" % (ma, mp, len(pairs))
        else:
            ma, mp, mtxt = np.nan, np.nan, "too few pairs (%d)" % len(pairs)
        rc = float(np.mean([D[c][d][2] > 0 for d in cc]))
        rk = float(np.mean([D[c][d][2] > 0 for d in ck]))
        flag = " *" if np.isfinite(mp) and mp < 0.01 else (" ." if np.isfinite(mp) and mp < 0.05 else "")
        print("%-9s %22s %26s %18s%s"
              % (c, "%.3f (p=%.2g)" % (pa, pp), mtxt, "%.0f%% / %.0f%%" % (100 * rc, 100 * rk), flag))
        rows.append({"cell": c, "set": args.set_name, "variant": args.variant,
                     "statistic": args.statistic, "contrast": "%s_vs_%s" % (args.case, args.control),
                     "n_case": na, "n_control": nb, "pooled_auroc": pa, "pooled_p": pp,
                     "matched_auroc": ma, "matched_p": mp, "n_pairs": len(pairs),
                     "responsive_case": rc, "responsive_control": rk})
    print()
    print("* matched p < 0.01, . matched p < 0.05. The matched column is the one to read:")
    print("the responsive share differs between cell lines, which is the confounder of report 08, 3.")
    if args.out:
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        with open(args.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
            w.writeheader()
            w.writerows(rows)
        print("\nwrote %s" % args.out, file=sys.stderr)


if __name__ == "__main__":
    main()
