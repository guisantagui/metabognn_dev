#!/usr/bin/env python3
"""Re-test a phenotype association after matching cases and controls on a covariate.

Generalises `matched_check.py`, which matches only on `log1p(n_sig)` responsiveness.
Here any numeric column of the phenotype table can be the matching covariate, which
matters for the Naia ATP screen: it carries a measured per-compound cytotoxicity
z-score, so viability can be controlled directly instead of through the `n_sig`
responsiveness proxy used in reports 07-09.

Why it is needed for that dataset specifically. ATP depletion and cell death are
confounded in principle, and glycolytic gene expression is a generic stress
response, so a glycolysis association could be "these compounds kill cells" rather
than anything metabolic. In this screen the two readouts are nearly independent
over all 1,199 compounds (Spearman = -0.083), but the mapped subset is small and
the ATP-down group is more cytotoxic than average, so it has to be checked on the
drugs actually scored rather than assumed from the full screen.

Matching is 1:1 greedy nearest-neighbour without replacement, hardest case first,
within a caliper, the same procedure as `matched_check.py`. Reports the pooled test,
each single-covariate match, and the joint match on both covariates (standardised
Euclidean distance), because responsiveness and cytotoxicity are different
confounders and controlling one does not control the other.
"""
import argparse, csv, os, sys

import numpy as np

try:
    from scipy import stats as sps
except Exception:                                   # noqa: BLE001
    sps = None


def read_phenotype(path):
    return {r['drug']: r for r in csv.DictReader(open(path), delimiter='\t')}


def read_scores(path, set_name, variant, statistic):
    out = {}
    for r in csv.DictReader(open(path), delimiter='\t'):
        if r['set'] != set_name or r['variant'] != variant:
            continue
        try:
            v = float(r[statistic])
        except (ValueError, KeyError):
            continue
        try:
            ns = int(r['max_n_sig'])
        except (ValueError, KeyError):
            ns = 0
        out[r['drug']] = {'group': r['group'], 'v': v, 'n_sig': ns}
    return out


def mw(a, b):
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    if len(a) < 3 or len(b) < 3 or sps is None:
        return np.nan, np.nan
    u, p = sps.mannwhitneyu(a, b, alternative='two-sided')
    return u / (len(a) * len(b)), p


def greedy(case, ctrl, caliper):
    """case/ctrl: list of (drug, coord-vector). Euclidean distance on standardised coords."""
    pool = list(ctrl)
    if not pool or not case:
        return []
    allc = np.array([c[1] for c in case] + [c[1] for c in pool], float)
    sd = allc.std(axis=0)
    sd[sd == 0] = 1.0
    cen = allc.mean(axis=0)
    norm = lambda v: (np.asarray(v, float) - cen) / sd
    centre = norm(np.median([c[1] for c in pool], axis=0))
    order = sorted(range(len(case)), key=lambda k: -np.linalg.norm(norm(case[k][1]) - centre))
    pairs = []
    for k in order:
        if not pool:
            break
        t = norm(case[k][1])
        d = [float(np.linalg.norm(norm(c[1]) - t)) for c in pool]
        j = int(np.argmin(d))
        if d[j] <= caliper:
            pairs.append((case[k][0], pool[j][0]))
            pool.pop(j)
    return pairs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--scores', required=True, help='per_drug_scores.tsv from mc_phenotype_assoc.py')
    ap.add_argument('--phenotype', required=True)
    ap.add_argument('--set', dest='set_name', required=True)
    ap.add_argument('--variant', default='all')
    ap.add_argument('--statistic', default='z')
    ap.add_argument('--case', required=True)
    ap.add_argument('--control', default='inactive')
    ap.add_argument('--covariate', action='append', default=[],
                    help='numeric phenotype column to match on; repeatable. '
                         'The literal name "n_sig" uses log1p(max_n_sig) from the scores file.')
    ap.add_argument('--caliper', type=float, default=0.5,
                    help='maximum standardised distance for a pair (default 0.5 SD)')
    ap.add_argument('--robustness', action='store_true',
                    help='also report leave-one-out, drop-the-k-most-extreme, the case listing '
                         'and the case MoA composition - the checks that matter when the case '
                         'group is small enough for a few drugs to carry the result')
    ap.add_argument('--drop-k', type=int, default=5,
                    help='with --robustness, how many most-extreme cases to drop (default 5)')
    args = ap.parse_args()

    pheno = read_phenotype(args.phenotype)
    sc = read_scores(args.scores, args.set_name, args.variant, args.statistic)
    if not sc:
        raise SystemExit('no rows for set=%s variant=%s statistic=%s in %s'
                         % (args.set_name, args.variant, args.statistic, args.scores))
    covs = args.covariate or ['n_sig']

    def coord(drug, names):
        out = []
        for c in names:
            if c == 'n_sig':
                out.append(float(np.log1p(sc[drug]['n_sig'])))
            else:
                try:
                    out.append(float(pheno[drug][c]))
                except (KeyError, ValueError, TypeError):
                    return None
        return out

    case = [d for d, r in sc.items() if r['group'] == args.case and np.isfinite(r['v'])]
    ctrl = [d for d, r in sc.items() if r['group'] == args.control and np.isfinite(r['v'])]
    print('%s / %s / %s   %s vs %s' % (args.set_name, args.variant, args.statistic,
                                       args.case, args.control))
    a, p = mw([sc[d]['v'] for d in case], [sc[d]['v'] for d in ctrl])
    print('  pooled                      AUROC %.3f  p=%.3g   n=%d/%d' % (a, p, len(case), len(ctrl)))

    for names in [[c] for c in covs] + ([covs] if len(covs) > 1 else []):
        cs = [(d, coord(d, names)) for d in case]
        ks = [(d, coord(d, names)) for d in ctrl]
        cs = [x for x in cs if x[1] is not None]
        ks = [x for x in ks if x[1] is not None]
        pairs = greedy(cs, ks, args.caliper)
        label = 'matched on ' + '+'.join(names)
        if len(pairs) < 10:
            print('  %-27s only %d pairs within the caliper - not reported' % (label, len(pairs)))
            continue
        a2, p2 = mw([sc[x]['v'] for x, _ in pairs], [sc[y]['v'] for _, y in pairs])
        # covariate balance before and after
        bal = []
        for i, c in enumerate(names):
            pre = np.mean([x[1][i] for x in cs]) - np.mean([x[1][i] for x in ks])
            post = (np.mean([dict(cs)[x][i] for x, _ in pairs])
                    - np.mean([dict(ks)[y][i] for _, y in pairs]))
            bal.append('%s %+.2f->%+.2f' % (c, pre, post))
        print('  %-27s AUROC %.3f  p=%.3g   %d pairs   [%s]'
              % (label, a2, p2, len(pairs), '; '.join(bal)))

    if args.robustness:
        import collections
        cv = [sc[d]['v'] for d in ctrl]
        print()
        print('  robustness of the pooled test')
        loo = []
        for d in case:
            a3, p3 = mw([sc[x]['v'] for x in case if x != d], cv)
            loo.append((p3, a3, d))
        loo.sort(reverse=True)
        print('    leave-one-out, worst 3 (p never exceeds the first):')
        for p3, a3, d in loo[:3]:
            print('      without %-26s AUROC %.3f  p=%.3g' % (d, a3, p3))
        k = min(args.drop_k, max(0, len(case) - 5))
        if k:
            kept = sorted(case, key=lambda x: sc[x]['v'])[k:]
            a4, p4 = mw([sc[x]['v'] for x in kept], cv)
            print('    dropping the %d most extreme cases: AUROC %.3f  p=%.3g  (n=%d)'
                  % (k, a4, p4, len(kept)))
        moa = collections.Counter((pheno.get(d, {}).get('moa') or '(none)') for d in case)
        print('    case MoA composition (top 6 of %d distinct): %s'
              % (len(moa), '; '.join('%s x%d' % (str(m)[:30], n) for m, n in moa.most_common(6))))
        print()
        print('  the %d case drugs by %s (lower = stronger shift):' % (len(case), args.statistic))
        for d in sorted(case, key=lambda x: sc[x]['v']):
            extra = ' '.join('%s=%s' % (c, pheno.get(d, {}).get(c, '?')) for c in covs if c != 'n_sig')
            print('      %-28s %7.2f  %s' % (d, sc[d]['v'], extra))


if __name__ == '__main__':
    main()
