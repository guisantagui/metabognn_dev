#!/usr/bin/env python3
"""Responsiveness-MATCHED version of the phenotype contrast.

`stratified_check.py` splits the drugs into responsive / non-responsive and
repeats the test inside each stratum. That controls the confounder but spends
most of the power (the responsive stratum of the MetChange run is 45 vs 44), and
it leaves the two strata to be reconciled by eye. This script instead removes the
confounder while keeping one estimate, two ways:

  matched       1:1 greedy nearest-neighbour matching of control drugs to case
                drugs on log1p(n_sig), without replacement, within --caliper.
                The matched control group then has, by construction, the same
                responsiveness distribution as the cases, and the plain
                Mann-Whitney on the matched pairs is unconfounded by it.

  residual      least-squares regression of the per-drug statistic on
                log1p(n_sig) over ALL drugs, then the same test on the residuals.
                Uses every drug instead of discarding unmatched controls, at the
                price of assuming the confounder acts linearly on that scale.

They fail in different ways - matching discards data, residualisation assumes a
form - so agreement between them is the thing to look for, not either alone.

Reads per_drug_scores.tsv, so no access to the label matrix is needed. n_sig must
be present in it (it is whenever the run had fit_info_nsig.csv or read a q-value
matrix; see --no-nsig in mc_phenotype_assoc.py).
"""
import argparse, collections, csv, os, sys

import numpy as np

try:
    from scipy import stats as sps
except Exception:                                      # noqa: BLE001
    sps = None


def mw(case, ctrl):
    case = np.asarray([x for x in case if np.isfinite(x)], float)
    ctrl = np.asarray([x for x in ctrl if np.isfinite(x)], float)
    if len(case) < 3 or len(ctrl) < 3 or sps is None:
        return np.nan, np.nan, len(case), len(ctrl)
    u, p = sps.mannwhitneyu(case, ctrl, alternative='two-sided')
    return u / (len(case) * len(ctrl)), p, len(case), len(ctrl)


def greedy_match(case_ns, ctrl_ns, caliper, seed=0):
    """1:1 greedy nearest-neighbour match on log1p(n_sig), no replacement.

    Cases are matched in order of how rare their responsiveness is (most extreme
    first), so the hardest ones get the pick of the controls rather than the
    leftovers.
    """
    rng = np.random.default_rng(seed)
    ctrl = [(i, v) for i, v in ctrl_ns]
    order = sorted(range(len(case_ns)), key=lambda k: -abs(
        case_ns[k][1] - np.median([v for _, v in ctrl_ns]) if ctrl_ns else 0.0))
    pairs = []
    for k in order:
        ci, cv = case_ns[k]
        if not ctrl:
            break
        d = [abs(v - cv) for _, v in ctrl]
        j = int(np.argmin(d))
        if d[j] <= caliper:
            pairs.append((ci, ctrl[j][0]))
            ctrl.pop(j)
    return pairs


def bh(ps):
    ps = np.asarray(ps, float)
    ok = np.isfinite(ps)
    q = np.full(len(ps), np.nan)
    p = ps[ok]
    if not len(p):
        return q
    o = p.argsort()
    r = p[o] * len(p) / np.arange(1, len(p) + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    out = np.empty(len(p))
    out[o] = np.clip(r, 0, 1)
    q[ok] = out
    return q


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--scores', required=True)
    ap.add_argument('--case', default='MMP_down')
    ap.add_argument('--control', default='inactive')
    ap.add_argument('--statistic', default='z')
    ap.add_argument('--caliper', type=float, default=0.25,
                    help='maximum |log1p(n_sig) difference| for a match')
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out')
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.scores), delimiter='\t'))
    if args.statistic not in rows[0]:
        raise SystemExit('%s has no column %r' % (args.scores, args.statistic))
    if 'max_n_sig' not in rows[0]:
        raise SystemExit('%s has no max_n_sig column' % args.scores)

    vals = collections.defaultdict(dict)
    meta = {}
    for r in rows:
        try:
            v = float(r[args.statistic])
        except ValueError:
            v = np.nan
        vals[(r['set'], r['variant'])][r['drug']] = v
        meta[r['drug']] = (r['group'], np.log1p(int(r['max_n_sig'])))

    drugs = sorted(meta)
    case_ns = [(d, meta[d][1]) for d in drugs if meta[d][0] == args.case]
    ctrl_ns = [(d, meta[d][1]) for d in drugs if meta[d][0] == args.control]
    if not case_ns or not ctrl_ns:
        raise SystemExit('one of the groups is empty')
    pairs = greedy_match(case_ns, ctrl_ns, args.caliper, args.seed)
    matched_case = [a for a, _ in pairs]
    matched_ctrl = [b for _, b in pairs]
    ns = lambda ds: np.array([meta[d][1] for d in ds])
    print('matching on log1p(n_sig), caliper %.2f: %d of %d cases matched'
          % (args.caliper, len(pairs), len(case_ns)))
    if pairs:
        print('  mean log1p(n_sig): cases %.3f, matched controls %.3f, all controls %.3f'
              % (ns(matched_case).mean(), ns(matched_ctrl).mean(),
                 ns([d for d, _ in ctrl_ns]).mean()))
        if sps is not None:
            print('  residual imbalance after matching: Mann-Whitney p = %.3g'
                  % sps.mannwhitneyu(ns(matched_case), ns(matched_ctrl),
                                     alternative='two-sided')[1])
    print()

    out = []
    hdr = '%-20s %-11s %18s %6s %6s %7s %10s' % (
        'set', 'variant', 'design', 'n_case', 'n_ctrl', 'AUROC', 'p')
    print(hdr)
    print('-' * len(hdr))
    for key in sorted(vals):
        v = vals[key]
        # 1. unmatched, for reference
        a, p, nc, nk = mw([v[d] for d, _ in case_ns], [v[d] for d, _ in ctrl_ns])
        rec = [('pooled (reference)', a, p, nc, nk)]
        # 2. matched
        a2, p2, nc2, nk2 = mw([v[d] for d in matched_case], [v[d] for d in matched_ctrl])
        rec.append(('responsiveness-matched', a2, p2, nc2, nk2))
        # 3. residual on log1p(n_sig)
        ds = [d for d in drugs if meta[d][0] in (args.case, args.control)
              and np.isfinite(v.get(d, np.nan))]
        if len(ds) > 10:
            y = np.array([v[d] for d in ds])
            x = np.array([meta[d][1] for d in ds])
            A = np.vstack([np.ones_like(x), x]).T
            beta, *_ = np.linalg.lstsq(A, y, rcond=None)
            res = y - A @ beta
            g = np.array([meta[d][0] for d in ds])
            a3, p3, nc3, nk3 = mw(res[g == args.case], res[g == args.control])
            rec.append(('residual on n_sig', a3, p3, nc3, nk3))
        for design, a_, p_, nc_, nk_ in rec:
            print('%-20s %-11s %18s %6d %6d %7s %10s'
                  % (key[0], key[1], design, nc_, nk_,
                     '%.3f' % a_ if np.isfinite(a_) else '-',
                     '%.3g' % p_ if np.isfinite(p_) else '-'))
            out.append({'set': key[0], 'variant': key[1], 'design': design,
                        'statistic': args.statistic, 'n_case': nc_, 'n_control': nk_,
                        'auroc': a_, 'p': p_})
        print()

    for design in {r['design'] for r in out}:
        sel = [r for r in out if r['design'] == design]
        for r, q in zip(sel, bh([r['p'] for r in sel])):
            r['bh_q'] = q

    if args.out:
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        with open(args.out, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=['set', 'variant', 'design', 'statistic',
                                               'n_case', 'n_control', 'auroc', 'p', 'bh_q'],
                               delimiter='\t')
            w.writeheader()
            w.writerows(out)
        print('wrote %s' % args.out, file=sys.stderr)


if __name__ == '__main__':
    main()
