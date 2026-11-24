#!/usr/bin/env python3
"""The within-drug, between-cell-line differential test on PRISM viability.

**This is the actual test; the per-cell group contrast in each run's
`set_results.tsv` is not** (report 10, 1.1). PRISM's phenotype is viability, and
viability is the confounder that dominated reports 07-09 - report 09, 9.6 turns on
it, the expression association dying under a cytotoxicity match while MetChange's
survived. You cannot match on cytotoxicity when cytotoxicity is the phenotype, so a
direct "do sensitive drugs score differently" test would likely light up every set,
glycan negative controls included.

THE DESIGN. For each drug d and ordered cell pair (X, Y):

    phenotype difference   dP = AUC(d, X) - AUC(d, Y)
    predictor difference   dS = setscore(d | X labels) - setscore(d | Y labels)

and correlate dP against dS across drugs. **Differencing removes the drug's general
potency by construction**: a compound that is simply toxic everywhere has dP ~ 0 and
contributes nothing to either side, so the confounder is cancelled rather than
matched away. What remains is strictly the selectivity question - when a drug hits
one line harder than another, do that line's MetChange labels say so?

PRE-REGISTERED DIRECTION (report 10, 1.2, written before the test was run). Lower
predicted production capacity for biomass precursors in X than in Y should mean MORE
growth inhibition in X, hence LOWER AUC in X. Both differences are negative
together, so **the expected correlation is POSITIVE**. A negative correlation is not
a result with the sign flipped; it contradicts the prediction.

CONTROLS, because a correlation alone is not enough:
  perm_p        drug identity of the dS vector is permuted, which breaks the
                drug-to-drug correspondence while preserving both marginal
                distributions. This is the null that matters.
  naive_r       the undifferenced correlation, setscore(X) against AUC(X), shown
                only as a foil: it is the confounded quantity the design exists to
                avoid, and it is expected to be larger and meaningless.
  selective     the same differential restricted to drugs whose |dP| exceeds
                --selective-cut. Pre-specified: a drug with no between-cell spread
                carries no signal for this design whichever way the prediction goes.
  bh_q          Benjamini-Hochberg over the whole set x variant x statistic x
                cell-pair family, which is large - with 5 cells there are 10
                unordered pairs and ~160 set x variant x statistic combinations.

Inputs are small: the phenotype table (a few hundred kB) and one
`per_drug_scores.tsv` per cell line (about 3 MB each). The GB-scale label matrices
are never touched here, so this runs equally well on the cluster or on a laptop
after downloading only those files.
"""
import argparse, csv, itertools, os, sys

import numpy as np

try:
    from scipy import stats as sps
except Exception:                                   # noqa: BLE001
    sps = None


def read_phenotype(path, metric):
    rows = list(csv.DictReader(open(path), delimiter='\t'))
    if not rows:
        raise SystemExit('empty phenotype table %s' % path)
    pref = '%s_' % metric
    cells = [c[len(pref):] for c in rows[0] if c.startswith(pref)]
    if not cells:
        raise SystemExit('no %s<CELL> columns in %s - was it built with --metric %s?'
                         % (pref, path, metric))
    out = {}
    for r in rows:
        d = {}
        for c in cells:
            try:
                d[c] = float(r['%s%s' % (pref, c)])
            except (ValueError, KeyError):
                pass
        if d:
            out[r['drug']] = d
    return out, sorted(cells)


def read_scores(path):
    """(set, variant, statistic) -> {drug: value}"""
    out = {}
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter='\t')
        stats = [c for c in rd.fieldnames if c not in
                 ('drug', 'group', 'subgroup', 'set', 'variant',
                  'n_perturbations', 'max_n_sig')]
        for r in rd:
            for st in stats:
                try:
                    v = float(r[st])
                except (ValueError, KeyError):
                    continue
                out.setdefault((r['set'], r['variant'], st), {})[r['drug']] = v
    return out


def bh(ps):
    ps = np.asarray(ps, float)
    ok = np.isfinite(ps)
    q = np.full(ps.shape, np.nan)
    p = ps[ok]
    if not len(p):
        return q
    order = np.argsort(p)
    m = len(p)
    adj = p[order] * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.clip(adj, 0, 1)
    q[ok] = out
    return q


def main():
    repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--phenotype', default=os.path.join(
        repo, 'results/mc_phenotype_assoc/phenotypes_ALL/prism_viability.tsv'))
    ap.add_argument('--scores-root', default=os.path.join(repo, 'results/mc_phenotype_assoc'),
                    help='holds <CELL>/<run-name>/per_drug_scores.tsv')
    ap.add_argument('--run-name', default='prism_viability')
    ap.add_argument('--metric', default='auc')
    ap.add_argument('--cells', default='', help='default: every cell with both a score file and a phenotype column')
    ap.add_argument('--statistic', default='z')
    ap.add_argument('--method', default='spearman', choices=['spearman', 'pearson'])
    ap.add_argument('--selective-cut', type=float, default=0.1,
                    help='|dP| above which a drug counts as cell-selective (default 0.1 AUC)')
    ap.add_argument('--min-drugs', type=int, default=30)
    ap.add_argument('--n-perm', type=int, default=2000)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--top', type=int, default=25, help='rows to print, ranked by BH q')
    ap.add_argument('--out')
    args = ap.parse_args()
    if sps is None:
        raise SystemExit('scipy is required')

    pheno, pheno_cells = read_phenotype(args.phenotype, args.metric)
    scores, missing = {}, []
    for c in (args.cells.replace(',', ' ').split() or pheno_cells):
        p = os.path.join(args.scores_root, c, args.run_name, 'per_drug_scores.tsv')
        if os.path.isfile(p):
            scores[c] = read_scores(p)
        else:
            missing.append(c)
    print('phenotype: %s  (%d drugs, cells %s)' % (args.phenotype, len(pheno), ', '.join(pheno_cells)))
    print('scored cells found: %s' % (', '.join(sorted(scores)) or 'NONE'))
    if missing:
        print('no per_drug_scores.tsv for: %s   (run launch_prism_slurm.sh for those)'
              % ', '.join(missing))
    if len(scores) < 2:
        raise SystemExit('\nthe differential needs at least TWO cell lines; '
                         'score another with: cell=<CELL> bash scripts/mc_phenotype_assoc/launch_prism_cell.sh')

    keys = sorted(set.intersection(*[{k for k in s if k[2] == args.statistic} for s in scores.values()]))
    pairs = list(itertools.combinations(sorted(scores), 2))
    print('%d set x variant combinations x %d cell pairs, statistic=%s, %s correlation'
          % (len(keys), len(pairs), args.statistic, args.method))
    corr = sps.spearmanr if args.method == 'spearman' else sps.pearsonr
    rng = np.random.default_rng(args.seed)

    rows = []
    for X, Y in pairs:
        for key in keys:
            sx, sy = scores[X][key], scores[Y][key]
            drugs = [d for d in pheno
                     if d in sx and d in sy and X in pheno[d] and Y in pheno[d]
                     and np.isfinite(sx[d]) and np.isfinite(sy[d])]
            if len(drugs) < args.min_drugs:
                continue
            dP = np.array([pheno[d][X] - pheno[d][Y] for d in drugs])
            dS = np.array([sx[d] - sy[d] for d in drugs])
            r, p = corr(dS, dP)
            # permutation null: break the drug correspondence, keep both margins
            null = np.array([corr(rng.permutation(dS), dP)[0] for _ in range(args.n_perm)])
            perm_p = (1 + np.sum(np.abs(null) >= abs(r))) / (1 + args.n_perm)
            sel = np.abs(dP) > args.selective_cut
            if sel.sum() >= args.min_drugs:
                r_sel, p_sel = corr(dS[sel], dP[sel])
            else:
                r_sel, p_sel = np.nan, np.nan
            # the confounded quantity, as a foil only
            naive = corr(np.array([sx[d] for d in drugs]),
                         np.array([pheno[d][X] for d in drugs]))[0]
            rows.append({'cell_x': X, 'cell_y': Y, 'set': key[0], 'variant': key[1],
                         'statistic': key[2], 'n_drugs': len(drugs), 'r_diff': r, 'p_diff': p,
                         'perm_p': perm_p, 'n_selective': int(sel.sum()),
                         'r_selective': r_sel, 'p_selective': p_sel, 'naive_r_in_x': naive})
    if not rows:
        raise SystemExit('no cell pair had >= %d shared drugs' % args.min_drugs)
    q = bh([r['p_diff'] for r in rows])
    for r, qq in zip(rows, q):
        r['bh_q'] = qq

    rows.sort(key=lambda r: (r['bh_q'], r['p_diff']))
    print()
    print('PRE-REGISTERED PREDICTION: r_diff POSITIVE (lower predicted capacity in X -> lower AUC in X).')
    print('A negative r_diff contradicts the prediction rather than supporting it with a flipped sign.')
    print()
    hdr = ('%-9s %-20s %-11s %6s %8s %8s %8s %8s %9s %9s'
           % ('pair', 'set', 'variant', 'n', 'r_diff', 'perm_p', 'BH q', 'r_sel', 'n_sel', 'naive_r'))
    print(hdr)
    print('-' * len(hdr))
    for r in rows[:args.top]:
        print('%-9s %-20s %-11s %6d %+8.3f %8.3g %8.3g %+8.3f %9d %+9.3f'
              % ('%s/%s' % (r['cell_x'], r['cell_y']), r['set'][:20], r['variant'],
                 r['n_drugs'], r['r_diff'], r['perm_p'], r['bh_q'],
                 r['r_selective'] if np.isfinite(r['r_selective']) else 0.0,
                 r['n_selective'], r['naive_r_in_x']))
    pos = sum(1 for r in rows if r['r_diff'] > 0)
    sig = [r for r in rows if r['bh_q'] < 0.05]
    print()
    print('%d of %d combinations have the predicted POSITIVE sign (%.0f%%; 50%% = no signal)'
          % (pos, len(rows), 100 * pos / len(rows)))
    print('%d pass BH q < 0.05, of which %d have the predicted positive sign'
          % (len(sig), sum(1 for r in sig if r['r_diff'] > 0)))
    print('median |naive_r| %.3f vs median |r_diff| %.3f - the naive value is the confounded one'
          % (np.median([abs(r['naive_r_in_x']) for r in rows]),
             np.median([abs(r['r_diff']) for r in rows])))
    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)) or '.', exist_ok=True)
        with open(args.out, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
            w.writeheader()
            w.writerows(rows)
        print('\nwrote %s' % args.out, file=sys.stderr)


if __name__ == '__main__':
    main()
