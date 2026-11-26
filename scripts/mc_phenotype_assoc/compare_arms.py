#!/usr/bin/env python3
"""Four-way table: MetChange vs producer-average vs global-residual vs expression.

Report 11, section 11. Reads the `set_results.tsv` that `mc_phenotype_assoc.py`
wrote for each arm and lines them up on the same (contrast, set, variant,
statistic) key, so the arms are compared on identical designs rather than on each
arm's own best row - which is the comparison that forking paths would otherwise
make for us.

RANDOM-SET p, NOT PARAMETRIC p. With MetChange features correlated at r=0.97
(section 8.1) the parametric Mann-Whitney null is not the relevant one: it assumes
features contribute independent information. `random_set_p` compares the set
against size-matched random sets drawn from the same correlated pool and is the
number to read. Its floor is 1/(n_random_sets+1), printed in the header, and a
value AT the floor means "not resolved by this many draws", not "this exact p".

The expression arm is scored on GENE sets, the other three on METABOLITE sets, so
set names match by construction but the underlying features differ. That is the
intended comparison: the same pathway, summarised from different representations.

VARIANTS DO NOT MATCH ACROSS ARMS. The metabolite arms carry compartment variants
(`all`, `cyto`, `mito`, and their `_nohub` forms); the gene arm carries `gene` and
`gene_nohub`, because a gene has no compartment. Keying on (set, variant) therefore
intersects to nothing. `--gene-arm` names the arms that should be matched on the set
name alone, taking their `gene` variant (or `gene_nohub` with --nohub), which is
unambiguous because there is exactly one of each per set. Only sets present in EVERY
arm are printed, so the table cannot silently become "each arm's best row".
"""

import argparse
import collections
import csv
import os
import sys

ARMS = [
    ("MetChange", "attene_mc_ref"),
    ("producer-avg", "attene_prodavg"),
    ("residual", "attene_resid"),
    ("expression", "attene_expr_ref"),
]
KEY = ("contrast", "set", "variant", "statistic")


def read_arm(path):
    """(contrast, set, variant, statistic) -> row."""
    if not os.path.exists(path):
        return None
    with open(path) as fh:
        return {tuple(r[k] for k in KEY): r for r in csv.DictReader(fh, delimiter="\t")}


def load(base, arms=ARMS):
    out, missing = {}, []
    for label, sub in arms:
        rows = read_arm(os.path.join(base, sub, "set_results.tsv"))
        if rows is None:
            missing.append((label, os.path.join(base, sub, "set_results.tsv")))
        else:
            out[label] = rows
    return out, missing


def table(data, statistic="z", contrast=None, sets=None, max_rows=40,
          order_by="MetChange", gene_arms=("expression",), nohub=False):
    """Print one row per (set, variant), one column per arm.

    Arms named in `gene_arms` are matched on the set name alone; see the module
    docstring on why variants cannot be keyed across arms.
    """
    labels = list(data)
    gene_variant = "gene_nohub" if nohub else "gene"

    def lookup(lab, key):
        """Row for this arm at this (contrast, set, variant, statistic) key."""
        if lab in gene_arms:
            return data[lab].get((key[0], key[1], gene_variant, key[3]))
        return data[lab].get(key)

    keys = set()
    for lab, rows in data.items():
        if lab in gene_arms:
            continue
        keys |= {k for k in rows if k[3] == statistic}
    if contrast:
        keys = {k for k in keys if k[0] == contrast}
    if sets:
        want = set(sets)
        keys = {k for k in keys if k[1] in want}
    common = [k for k in keys if all(lookup(lab, k) is not None for lab in labels)]
    if not common:
        print("  no set is present in all %d arms for statistic=%s" % (len(labels), statistic))
        partial = collections.Counter()
        for k in keys:
            for lab in labels:
                if lookup(lab, k) is not None:
                    partial[lab] += 1
        print("  rows per arm at this statistic: %s" % dict(partial))
        return []

    def sortkey(k):
        r = lookup(order_by if order_by in data else labels[0], k)
        return float(r["random_set_p"]) if r else 1.0

    common.sort(key=sortkey)
    floors = {lab: 1.0 / (int(next(iter(rows.values()))["n_random_sets"]) + 1)
              for lab, rows in data.items()}
    print("  statistic=%s   contrast=%s   (gene arms matched on set name, variant=%s)"
          % (statistic, contrast or "ALL", gene_variant))
    print("  random-set p floors: %s"
          % ", ".join("%s %.4g" % (k, v) for k, v in floors.items()))
    print()
    print("  %-20s %-11s" % ("set", "variant")
          + "".join(" | %-21s" % lab for lab in labels))
    print("  %-20s %-11s" % ("", "")
          + "".join(" | %-21s" % "AUROC   rand     n" for _ in labels))
    print("  " + "-" * (33 + 24 * len(labels)))
    for k in common[:max_rows]:
        line = "  %-20s %-11s" % (k[1][:20], k[2][:11])
        for lab in labels:
            r = lookup(lab, k)
            line += " | %6.3f %7.4g %5s" % (
                float(r["group_auroc"]), float(r["random_set_p"]), r["n_case"])
        print(line)
    return common


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--base", default="results/mc_phenotype_assoc/HEPG2")
    ap.add_argument("--statistic", default="z", choices=["z", "auroc"])
    ap.add_argument("--contrast", default=None,
                    help="e.g. MMP_down_vs_inactive; default prints each in turn")
    ap.add_argument("--sets", nargs="*", default=None)
    ap.add_argument("--max-rows", type=int, default=40)
    ap.add_argument("--gene-arm", nargs="*", default=["expression"],
                    help="arms matched on set name alone (gene sets have no compartment)")
    ap.add_argument("--nohub", action="store_true",
                    help="use the gene_nohub variant for the gene arms")
    a = ap.parse_args()
    data, missing = load(a.base)
    for label, path in missing:
        print("MISSING arm %-14s %s" % (label, path), file=sys.stderr)
    if len(data) < 2:
        raise SystemExit("need at least two arms; run reproduce_report11.sh 9 10 first")
    contrasts = ([a.contrast] if a.contrast
                 else sorted({k[0] for rows in data.values() for k in rows}))
    for con in contrasts:
        print()
        table(data, statistic=a.statistic, contrast=con, sets=a.sets,
              max_rows=a.max_rows, gene_arms=tuple(a.gene_arm), nohub=a.nohub)


if __name__ == "__main__":
    main()
