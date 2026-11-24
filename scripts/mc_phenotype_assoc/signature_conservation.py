#!/usr/bin/env python3
"""How conserved is a drug's MetChange signature between two cell lines?

The question this answers. Reports 09 and 10 contain two results that look
contradictory until this is measured:

  - report 09, section 3: for a HepG2-measured membrane-potential phenotype, only
    HEPG2 labels predict it (AUROC 0.350); the other seven cell lines are null.
  - report 09, section 9 and report 10, section 5: HEPG2 labels nonetheless predict
    an SH-SY5Y ATP phenotype, and lean toward predicting AC16 membrane potential.

A compound's mechanism does not change between cell lines, and the TCA cycle,
glycolysis and OXPHOS are common to all of them, so **some transfer is the expected
baseline rather than a finding**. The real question is quantitative: how much of a
MetChange signature is compound-driven (hence shared across cell lines) and how much
is cell-context-driven? If signatures are strongly conserved, transfer is trivial and
it is section 3's cell-specificity that needs explaining. If they are weakly
conserved, section 3 is unsurprising and the transfer is the notable part.

What it computes, per drug present in both cell lines:

  same-drug r     Pearson/Spearman correlation between the drug's signature in cell
                  A and in cell B, over the shared metabolites.
  different-drug r the same, for mismatched drug pairs.

**The second is the control that matters.** A high same-drug correlation proves
nothing on its own: if some metabolites always move under any perturbation, every
pair of signatures correlates. Conservation means same-drug correlation EXCEEDS
different-drug correlation; the gap, not the level, is the quantity of interest.

Signatures are collapsed to one vector per drug by averaging that drug's
perturbations within each cell line, as elsewhere in this package, because doses and
times of one drug are not independent.

**RECORD THE LABEL MATRIX mtime WITH ANY RESULT.** This script reads
`results/labels/<CELL>/...` directly, and MetChange labels exist in more than one
generation distinguished only by date: as of 2026-11-24 the **2026-11-12 files are
the corrected ones** and a superseded buggy generation was dated 2026-11-02. Nothing
in a file's name or content says which it is, and reading a superseded matrix gives a
confident wrong number with no error. Check the mtime and state it in the write-up.
(The first run of this script, 2026-11-24 HEPG2 vs A549, was briefly withdrawn on a
mistaken provenance worry and then reinstated - it had used the correct Nov-12
files.)
"""
import argparse, csv, os, sys

import numpy as np
import pandas as pd


def read_fit(path):
    """perturbation -> drug."""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            out[r['perturbation']] = r.get('drug', '')
    return out


def drug_matrix(label_dir, drugs, matrix, fit_info, max_perts, rows, seed):
    fit = read_fit(os.path.join(label_dir, fit_info))
    want = {}
    for pert, drug in fit.items():
        if drug in drugs:
            want.setdefault(drug, []).append(pert)
    rng = np.random.default_rng(seed)
    cols = []
    for drug, perts in want.items():
        if len(perts) > max_perts:
            perts = list(rng.choice(perts, max_perts, replace=False))
            want[drug] = perts
        cols += perts
    path = os.path.join(label_dir, matrix)
    header = pd.read_csv(path, nrows=0).columns.tolist()
    idcol = header[0]
    present = [c for c in cols if c in set(header)]
    print('  %s: %d drugs, %d perturbation columns read' % (label_dir, len(want), len(present)),
          file=sys.stderr)
    df = pd.read_csv(path, usecols=[idcol] + present, index_col=0)
    if rows and rows < len(df):
        df = df.iloc[np.sort(rng.choice(len(df), rows, replace=False))]
    out = pd.DataFrame({d: df[[p for p in ps if p in df.columns]].mean(axis=1)
                        for d, ps in want.items() if any(p in df.columns for p in ps)})
    return out


def main():
    repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--cell-a', default='HEPG2')
    ap.add_argument('--cell-b', default='A549')
    ap.add_argument('--labels-root', default=os.path.join(repo, 'results/labels'))
    ap.add_argument('--sub-path', default='metchange/rank',
                    help="e.g. 'metchange/rank' or 'expr'")
    ap.add_argument('--matrix', default='log2FC.csv')
    ap.add_argument('--fit-info', default='fit_info_nsig.csv')
    ap.add_argument('--n-drugs', type=int, default=250,
                    help='how many shared drugs to sample (these matrices are GB-scale)')
    ap.add_argument('--max-perts', type=int, default=3,
                    help='cap on perturbations averaged per drug, to bound the columns read')
    ap.add_argument('--rows', type=int, default=0, help='subsample metabolites (0 = all)')
    ap.add_argument('--n-pairs', type=int, default=20000, help='mismatched pairs for the control')
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out')
    args = ap.parse_args()

    dirs = {c: os.path.join(args.labels_root, c, args.sub_path) for c in (args.cell_a, args.cell_b)}
    for c, d in dirs.items():
        if not os.path.exists(os.path.join(d, args.matrix)):
            raise SystemExit('no %s in %s - is %s available locally?' % (args.matrix, d, c))

    fits = {c: read_fit(os.path.join(d, args.fit_info)) for c, d in dirs.items()}
    shared = sorted(set(fits[args.cell_a].values()) & set(fits[args.cell_b].values()) - {''})
    rng = np.random.default_rng(args.seed)
    if args.n_drugs and args.n_drugs < len(shared):
        shared = sorted(rng.choice(shared, args.n_drugs, replace=False))
    print('drugs shared by %s and %s: %d used' % (args.cell_a, args.cell_b, len(shared)),
          file=sys.stderr)

    mats = {c: drug_matrix(dirs[c], set(shared), args.matrix, args.fit_info,
                           args.max_perts, args.rows, args.seed) for c in dirs}
    feats = mats[args.cell_a].index.intersection(mats[args.cell_b].index)
    drugs = [d for d in mats[args.cell_a].columns if d in set(mats[args.cell_b].columns)]
    A = mats[args.cell_a].loc[feats, drugs].to_numpy(float)
    B = mats[args.cell_b].loc[feats, drugs].to_numpy(float)
    print('%d metabolites x %d drugs in both' % (len(feats), len(drugs)), file=sys.stderr)

    def z(M):
        M = M - np.nanmean(M, axis=0, keepdims=True)
        sd = np.nanstd(M, axis=0, keepdims=True)
        sd[sd == 0] = np.nan
        return M / sd

    Az, Bz = z(A), z(B)
    n = len(feats)
    same = np.nansum(Az * Bz, axis=0) / n
    i = rng.integers(0, len(drugs), args.n_pairs)
    j = rng.integers(0, len(drugs), args.n_pairs)
    keep = i != j
    diff = np.nansum(Az[:, i[keep]] * Bz[:, j[keep]], axis=0) / n

    def q(v):
        v = v[np.isfinite(v)]
        return np.median(v), np.percentile(v, 25), np.percentile(v, 75), v

    ms, s25, s75, sv = q(same)
    md, d25, d75, dv = q(diff)
    print()
    print('cross-cell signature correlation, %s vs %s (%s)'
          % (args.cell_a, args.cell_b, args.sub_path))
    print('  SAME drug      median r = %+.3f   IQR %+.3f .. %+.3f   (n=%d drugs)' % (ms, s25, s75, len(sv)))
    print('  DIFFERENT drug median r = %+.3f   IQR %+.3f .. %+.3f   (n=%d pairs)' % (md, d25, d75, len(dv)))
    print('  gap (same - different medians) = %+.3f' % (ms - md))
    try:
        from scipy import stats as sps
        u, p = sps.mannwhitneyu(sv, dv, alternative='two-sided')
        print('  same vs different: AUROC %.3f  p=%.3g' % (u / (len(sv) * len(dv)), p))
        print('  fraction of same-drug r above the different-drug 95th percentile: %.1f%%'
              % (100 * np.mean(sv > np.percentile(dv, 95))))
    except Exception:                                   # noqa: BLE001
        pass
    print()
    print('  Read the GAP, not the level: if some metabolites move under any perturbation,')
    print('  every pair of signatures correlates. Conservation is same exceeding different.')
    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)) or '.', exist_ok=True)
        pd.DataFrame({'drug': drugs, 'same_drug_r': same}).to_csv(args.out, sep='\t', index=False)
        print('\nwrote %s' % args.out, file=sys.stderr)


if __name__ == '__main__':
    main()
