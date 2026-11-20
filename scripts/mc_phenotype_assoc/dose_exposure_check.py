#!/usr/bin/env python3
"""Did LINCS actually profile the drug at a concentration where the phenotype occurs?

The phenotype call and the transcriptional profile come from different
experiments at different concentrations. Tox21 qHTS titrates to ~92 uM and calls
a compound active if its curve fits anywhere in that range; LINCS profiled these
compounds at fixed doses (39% of the relevant HepG2 perturbations at 10 uM). So a
compound with AC50 = 30 uM is "MMP-active" while its LINCS signature was taken
well below the concentration at which the membrane potential actually collapses.
Measured here: 48% of the MMP-down cases were profiled below their own AC50.

That makes the case group diluted by construction, and gives a falsifiable
prediction. If the TCA association is a real metabolic response to the phenotype:

  - cases profiled AT OR ABOVE their AC50 should show a STRONGER effect than the
    pooled estimate;
  - cases profiled BELOW their AC50 should look like inactives.

If instead both halves show it equally, that favours a compound-property
explanation, because a chemotype does not care at what dose the drug was
profiled.

Exposure is assigned per drug from its HIGHEST profiled dose: did any of its
profiles reach the phenotype concentration. 71% of these drugs have exactly one
perturbation, so for most the assignment is exact; for the 28% profiled at 4-7
doses the averaged score mixes exposed and unexposed perturbations, which dilutes
the contrast. A per-perturbation version would need the per-column statistics in
set_member_values.npz rather than per_drug_scores.tsv.

**Two limitations that bound what this test can show**, both of which weaken it
relative to how it was first proposed:

  1. The LINCS dose does not really vary - the interquartile spread of the maximum
     profiled dose across these drugs is 0.00 decades (essentially all 10 uM),
     while AC50 spans 0.86 decades. So log10(dose / AC50) is driven almost
     entirely by AC50, and this is therefore NOT independent of the potency
     correlation: it is close to the same test dichotomised at 10 uM.
  2. The MMP AC50 comes from a 1-5 h ratiometric dye assay. It is not a threshold
     for a 6-24 h transcriptional or metabolic change, and a drug may perturb
     metabolism at concentrations below full potential collapse. So "profiled
     below its AC50" does not mean "no metabolic effect was there to detect", and
     a failed prediction is correspondingly weaker evidence.

Note the direction convention: the group effect has AUROC < 0.5 (cases score
lower), so a STRONGER effect means a SMALLER AUROC.
"""
import argparse, collections, csv, os, sys

import numpy as np

try:
    from scipy import stats as sps
except Exception:                                      # noqa: BLE001
    sps = None


def mw(a, b):
    a = np.asarray([x for x in a if np.isfinite(x)], float)
    b = np.asarray([x for x in b if np.isfinite(x)], float)
    if len(a) < 3 or len(b) < 3 or sps is None:
        return np.nan, np.nan, len(a), len(b)
    u, p = sps.mannwhitneyu(a, b, alternative='two-sided')
    return u / (len(a) * len(b)), p, len(a), len(b)


def boot_ci(a, b, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    a = np.asarray([x for x in a if np.isfinite(x)], float)
    b = np.asarray([x for x in b if np.isfinite(x)], float)
    if len(a) < 3 or len(b) < 3:
        return np.nan, np.nan
    v = []
    for _ in range(n):
        x = rng.choice(a, len(a), replace=True)
        y = rng.choice(b, len(b), replace=True)
        v.append(mw(x, y)[0])
    return float(np.nanpercentile(v, 2.5)), float(np.nanpercentile(v, 97.5))


def main():
    repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--scores', required=True)
    ap.add_argument('--phenotype', required=True)
    ap.add_argument('--labels', required=True, help='for fit_info_nsig.csv (doses)')
    ap.add_argument('--set', default='tca')
    ap.add_argument('--variant', default='mito_nohub')
    ap.add_argument('--statistic', default='z')
    ap.add_argument('--case', default='MMP_down')
    ap.add_argument('--control', default='inactive')
    ap.add_argument('--fold', type=float, default=1.0,
                    help='a perturbation counts as exposed if dose >= fold x AC50')
    ap.add_argument('--out')
    args = ap.parse_args()

    pheno = {r['drug']: r for r in csv.DictReader(open(args.phenotype), delimiter='\t')}
    doses = collections.defaultdict(list)
    for r in csv.DictReader(open(os.path.join(args.labels, 'fit_info_nsig.csv'))):
        try:
            doses[r['drug']].append((r['perturbation'], float(r['dose'])))
        except (ValueError, KeyError):
            pass

    per_drug, per_pert = {}, collections.defaultdict(list)
    for r in csv.DictReader(open(args.scores), delimiter='\t'):
        if r['set'] != args.set or r['variant'] != args.variant:
            continue
        try:
            v = float(r[args.statistic])
        except (ValueError, KeyError):
            v = np.nan
        per_drug[r['drug']] = (r['group'], v)
    if not per_drug:
        raise SystemExit('set %s/%s not in %s' % (args.set, args.variant, args.scores))

    def ac50(d):
        s = (pheno.get(d, {}).get('potency_uM') or '').strip()
        try:
            return float(s)
        except ValueError:
            return None

    ctrl = [v for d, (g, v) in per_drug.items() if g == args.control and np.isfinite(v)]
    cases = [(d, v) for d, (g, v) in per_drug.items() if g == args.case and np.isfinite(v)]

    exposed, sub, noac = [], [], []
    detail = []
    for d, v in cases:
        a = ac50(d)
        dl = [x[1] for x in doses.get(d, [])]
        if a is None or not dl:
            noac.append(v)
            cls = 'no AC50 or no dose'
        elif max(dl) >= args.fold * a:
            exposed.append(v)
            cls = 'exposed'
        else:
            sub.append(v)
            cls = 'sub-threshold'
        detail.append({'drug': d, 'statistic_value': v, 'ac50_uM': a,
                       'max_lincs_dose_uM': max(dl) if dl else None, 'exposure_class': cls,
                       'log10_dose_over_ac50': (np.log10(max(dl) / a) if a and dl else None)})

    rows = []
    print('%s %s / %s / %s vs %s   (fold = %.3g)'
          % (args.set, args.variant, args.statistic, args.case, args.control, args.fold))
    print('controls: %d drugs' % len(ctrl))
    for lab, vals in (('ALL cases (pooled)', [v for _, v in cases]),
                      ('cases profiled AT/ABOVE AC50', exposed),
                      ('cases profiled BELOW AC50', sub),
                      ('cases without AC50/dose', noac)):
        a_, p_, n1, n2 = mw(vals, ctrl)
        lo, hi = boot_ci(vals, ctrl)
        print('  %-30s n=%4d  AUROC %s [%s, %s]  p = %s'
              % (lab, len(vals),
                 '%.3f' % a_ if np.isfinite(a_) else '-',
                 '%.3f' % lo if np.isfinite(lo) else '-',
                 '%.3f' % hi if np.isfinite(hi) else '-',
                 '%.3g' % p_ if np.isfinite(p_) else '-'))
        rows.append({'comparison': lab, 'n_case': len(vals), 'n_control': len(ctrl),
                     'auroc': a_, 'ci_lo': lo, 'ci_hi': hi, 'p': p_})

    # exposed vs sub-threshold directly: the prediction is exposed < sub-threshold
    a_, p_, n1, n2 = mw(exposed, sub)
    print('  %-30s n=%4d vs %-4d AUROC %s  p = %s'
          % ('exposed vs sub-threshold', n1, n2,
             '%.3f' % a_ if np.isfinite(a_) else '-',
             '%.3g' % p_ if np.isfinite(p_) else '-'))
    rows.append({'comparison': 'exposed vs sub-threshold (cases only)',
                 'n_case': n1, 'n_control': n2, 'auroc': a_, 'ci_lo': np.nan,
                 'ci_hi': np.nan, 'p': p_})

    # continuous version: within cases, score vs log10(dose / AC50)
    xs = [r['log10_dose_over_ac50'] for r in detail if r['log10_dose_over_ac50'] is not None]
    ys = [r['statistic_value'] for r in detail if r['log10_dose_over_ac50'] is not None]
    if len(xs) > 10 and sps is not None:
        rho, prho = sps.spearmanr(xs, ys)
        print('  within cases: Spearman(log10 dose/AC50, %s) = %+.3f, p = %.3g  (n=%d)'
              % (args.statistic, rho, prho, len(xs)))
        print('      negative rho = more exposed drugs score LOWER = phenotype-consistent')
        rows.append({'comparison': 'Spearman(log10 dose/AC50, statistic) within cases',
                     'n_case': len(xs), 'n_control': 0, 'auroc': rho,
                     'ci_lo': np.nan, 'ci_hi': np.nan, 'p': prho})

    if args.out:
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        with open(args.out, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=['comparison', 'n_case', 'n_control',
                                               'auroc', 'ci_lo', 'ci_hi', 'p'],
                               delimiter='\t')
            w.writeheader()
            w.writerows(rows)
        with open(args.out.replace('.tsv', '_perdrug.tsv'), 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=['drug', 'statistic_value', 'ac50_uM',
                                               'max_lincs_dose_uM', 'log10_dose_over_ac50',
                                               'exposure_class'], delimiter='\t')
            w.writeheader()
            w.writerows(detail)
        print('wrote %s' % args.out, file=sys.stderr)


if __name__ == '__main__':
    main()
