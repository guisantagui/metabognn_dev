#!/usr/bin/env python3
"""Is the phenotype association explained by MetChange responsiveness?

MMP-active compounds are bioactive compounds, and in HepG2 they are far more
likely to have a responsive MetChange signature at all (31% vs 9% of inactives
have a perturbation with a significant feature). Both per-drug statistics are
normalised inside each perturbation, so a larger overall response does not
mechanically raise a set's score - but if responsive signatures have a different
*shape* (set members systematically ranked low relative to other features), then
responsiveness alone would reproduce the association.

So the contrast is repeated inside each responsiveness stratum. An effect present
in both strata is not explained by responsiveness; an effect only in the pooled
data is confounded. Also reports the effect of responsiveness *within* the
control group, which is the confounder measured directly, and repeats the
contrast on 24 h perturbations only (the groups differ slightly in time
composition: 71% vs 62% at 24 h).

Reads per_drug_scores.tsv, so it needs no access to the label matrix.
"""
import argparse, collections, csv, os, sys

import numpy as np

try:
    from scipy import stats as sps
except Exception:                                      # noqa: BLE001
    sps = None


def auroc_p(case, ctrl):
    case = np.asarray([x for x in case if np.isfinite(x)], float)
    ctrl = np.asarray([x for x in ctrl if np.isfinite(x)], float)
    if len(case) < 3 or len(ctrl) < 3:
        return np.nan, np.nan, len(case), len(ctrl)
    if sps is None:
        return np.nan, np.nan, len(case), len(ctrl)
    u, p = sps.mannwhitneyu(case, ctrl, alternative='two-sided')
    return u / (len(case) * len(ctrl)), p, len(case), len(ctrl)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--scores', required=True, help='per_drug_scores.tsv')
    ap.add_argument('--case', default='MMP_down')
    ap.add_argument('--control', default='inactive')
    ap.add_argument('--statistic', default='auroc', choices=['auroc', 'z'])
    ap.add_argument('--sets', default='tca_endpoints:cyto_nohub,tca:mito_nohub,'
                                     'pyruvate:mito_nohub,glycolysis:cyto_nohub,'
                                     'redox_couples:all,energy_currency:mito,oxphos:mito',
                    help='comma-separated set:variant pairs to report')
    ap.add_argument('--out')
    args = ap.parse_args()

    rows = list(csv.DictReader(open(args.scores), delimiter='\t'))
    want = [tuple(s.split(':', 1)) for s in args.sets.split(',')]

    by = collections.defaultdict(dict)
    meta = {}
    for r in rows:
        try:
            v = float(r[args.statistic])
        except (ValueError, KeyError):
            v = np.nan
        by[(r['set'], r['variant'])][r['drug']] = v
        meta[r['drug']] = (r['group'], int(r['max_n_sig']))

    out = []
    header = ('%-20s %-11s %-22s %6s %6s %7s %10s' %
              ('set', 'variant', 'stratum', 'n_case', 'n_ctrl', 'AUROC', 'p'))
    print(header)
    print('-' * len(header))
    for key in want:
        vals = by.get(key)
        if not vals:
            print('%-20s %-11s  (not in the score table)' % key, file=sys.stderr)
            continue
        strata = {
            'pooled':          lambda ns: True,
            'responsive':      lambda ns: ns > 0,
            'non-responsive':  lambda ns: ns == 0,
        }
        for sname, keep in strata.items():
            case = [v for d, v in vals.items()
                    if meta[d][0] == args.case and keep(meta[d][1])]
            ctrl = [v for d, v in vals.items()
                    if meta[d][0] == args.control and keep(meta[d][1])]
            a, p, nc, nk = auroc_p(case, ctrl)
            print('%-20s %-11s %-22s %6d %6d %7s %10s'
                  % (key[0], key[1], sname, nc, nk,
                     '%.3f' % a if np.isfinite(a) else '-',
                     '%.3g' % p if np.isfinite(p) else '-'))
            out.append({'set': key[0], 'variant': key[1], 'stratum': sname,
                        'n_case': nc, 'n_control': nk, 'auroc': a, 'p': p,
                        'statistic': args.statistic})
        # the confounder measured directly: responsiveness inside the control group
        resp = [v for d, v in vals.items() if meta[d][0] == args.control and meta[d][1] > 0]
        nonr = [v for d, v in vals.items() if meta[d][0] == args.control and meta[d][1] == 0]
        a, p, nc, nk = auroc_p(resp, nonr)
        print('%-20s %-11s %-22s %6d %6d %7s %10s'
              % (key[0], key[1], 'resp vs non (controls)', nc, nk,
                 '%.3f' % a if np.isfinite(a) else '-',
                 '%.3g' % p if np.isfinite(p) else '-'))
        out.append({'set': key[0], 'variant': key[1], 'stratum': 'responsiveness_within_controls',
                    'n_case': nc, 'n_control': nk, 'auroc': a, 'p': p,
                    'statistic': args.statistic})
        print()

    if args.out:
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        with open(args.out, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=['set', 'variant', 'stratum', 'statistic',
                                               'n_case', 'n_control', 'auroc', 'p'],
                               delimiter='\t')
            w.writeheader()
            w.writerows(out)
        print('wrote %s' % args.out, file=sys.stderr)


if __name__ == '__main__':
    main()
