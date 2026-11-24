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
import argparse, csv, itertools, os, re, sys

import numpy as np

try:
    from scipy import stats as sps
except Exception:                                   # noqa: BLE001
    sps = None


def read_phenotype(path, metric, moa_out=None):
    rows = list(csv.DictReader(open(path), delimiter='\t'))
    if moa_out is not None:
        for r in rows:
            moa_out[r['drug']] = r.get('moa', '') or ''
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
    ap.add_argument('--moa-regex', default='',
                    help='case-insensitive pattern on the phenotype table\'s `moa` column. When given, every '
                         'combination is additionally reported for the MATCHING drugs and for the COMPLEMENT '
                         'separately. This matters more than any other option here: growth inhibition is '
                         'achieved by many mechanisms, and a hypothesis about metabolite production applies '
                         'only to drugs that act through metabolism. In the PRISM repurposing set just 45 of '
                         '821 drugs (5%) have a metabolic MoA - the rest are kinase inhibitors, tubulin '
                         'poisons, HDAC inhibitors, receptor ligands. Simulation: with a TRUE r of 0.8 among '
                         'those 45, the pooled 821-drug test reports r ~ 0.045 and reaches p<0.05 only 24% of '
                         'the time, so the pooled test cannot distinguish a strong mechanism-specific effect '
                         'from none at all. The complement is the built-in control: the hypothesis predicts '
                         'an effect in the matching group and NOT in the complement.')
    ap.add_argument('--moa-min-drugs', type=int, default=20,
                    help='minimum drugs in a MoA subgroup before it is reported (default 20)')
    ap.add_argument('--standardise-within-cell', action='store_true',
                    help='z-score each cell line\'s scores across drugs, per set, BEFORE differencing. '
                         'Differencing only cancels the confounder if neither cell dominates the '
                         'difference; when one cell\'s scores are more dispersed, dS leans towards it and '
                         'r_diff inherits that cell\'s naive within-cell correlation. Observed: '
                         'corr(r_diff, naive_r_in_x) = +0.40 overall and +0.45 for A375, which has the '
                         'largest dispersion of the five lines (SD 0.328 vs 0.263-0.303), and every '
                         'BH-significant result was a pair CONTAINING A375. Use this to test whether an '
                         'effect survives once that asymmetry is removed.')
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

    moa = {}
    pheno, pheno_cells = read_phenotype(args.phenotype, args.metric, moa)
    moa_re = re.compile(args.moa_regex, re.I) if args.moa_regex else None
    if moa_re is not None:
        n_match = sum(1 for d in pheno if moa_re.search(moa.get(d, '')))
        print('MoA filter %r matches %d of %d drugs (%.0f%%); the complement is the control'
              % (args.moa_regex, n_match, len(pheno), 100 * n_match / max(1, len(pheno))))
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
    if args.standardise_within_cell:
        for c, sc in scores.items():
            for key in list(sc):
                vals = np.array([v for v in sc[key].values() if np.isfinite(v)], float)
                if len(vals) < 3:
                    continue
                mu, sd = vals.mean(), vals.std()
                if sd == 0:
                    continue
                sc[key] = {d: (v - mu) / sd for d, v in sc[key].items()}
        print('scores z-scored within each cell line per set before differencing '
              '(--standardise-within-cell)')
    pairs = list(itertools.combinations(sorted(scores), 2))
    print('%d set x variant combinations x %d cell pairs, statistic=%s, %s correlation'
          % (len(keys), len(pairs), args.statistic, args.method))
    corr = sps.spearmanr if args.method == 'spearman' else sps.pearsonr
    rng = np.random.default_rng(args.seed)

    rows, degenerate = [], []
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
            # A set whose members are all absent from both label matrices gives a
            # constant (or all-nan) difference vector. Guard it: with a non-finite r,
            # `|null| >= |r|` is False for every permutation, so the permutation p-value
            # would come out as 1/(1+n_perm) - a spuriously SIGNIFICANT value for an
            # input that carries no information at all. Caught by the dS == 0 smoke test.
            if (np.nanstd(dS) == 0 or np.nanstd(dP) == 0
                    or not np.isfinite(dS).any() or not np.isfinite(dP).any()):
                degenerate.append((X, Y, key))
                continue
            r, p = corr(dS, dP)
            if not np.isfinite(r):
                degenerate.append((X, Y, key))
                continue
            # permutation null: break the drug correspondence, keep both margins
            null = np.array([corr(rng.permutation(dS), dP)[0] for _ in range(args.n_perm)])
            null = null[np.isfinite(null)]
            perm_p = ((1 + np.sum(np.abs(null) >= abs(r))) / (1 + len(null))
                      if len(null) else np.nan)
            sel = np.abs(dP) > args.selective_cut
            if sel.sum() >= args.min_drugs and np.nanstd(dS[sel]) > 0:
                r_sel, p_sel = corr(dS[sel], dP[sel])
            else:
                r_sel, p_sel = np.nan, np.nan
            # the confounded quantity, as a foil only
            naive = corr(np.array([sx[d] for d in drugs]),
                         np.array([pheno[d][X] for d in drugs]))[0]
            # MoA strata: the matching group is where the hypothesis applies at all,
            # the complement is the built-in control.
            r_moa = p_moa = r_oth = p_oth = np.nan
            n_moa = n_oth = 0
            if moa_re is not None:
                m = np.array([bool(moa_re.search(moa.get(d, ''))) for d in drugs])
                n_moa, n_oth = int(m.sum()), int((~m).sum())
                if n_moa >= args.moa_min_drugs and np.nanstd(dS[m]) > 0 and np.nanstd(dP[m]) > 0:
                    r_moa, p_moa = corr(dS[m], dP[m])
                if n_oth >= args.moa_min_drugs and np.nanstd(dS[~m]) > 0 and np.nanstd(dP[~m]) > 0:
                    r_oth, p_oth = corr(dS[~m], dP[~m])
            rows.append({'cell_x': X, 'cell_y': Y, 'set': key[0], 'variant': key[1],
                         'statistic': key[2], 'n_drugs': len(drugs), 'r_diff': r, 'p_diff': p,
                         'perm_p': perm_p, 'n_selective': int(sel.sum()),
                         'r_selective': r_sel, 'p_selective': p_sel, 'naive_r_in_x': naive,
                         'n_moa': n_moa, 'r_moa': r_moa, 'p_moa': p_moa,
                         'n_other': n_oth, 'r_other': r_oth, 'p_other': p_oth})
    if degenerate:
        print('\nSKIPPED %d combination(s) with a constant or all-missing difference vector '
              '(no information; reporting them would invent a p-value):' % len(degenerate))
        for X, Y, key in degenerate[:8]:
            print('   %s/%s  %s %s %s' % (X, Y, key[0], key[1], key[2]))
        if len(degenerate) > 8:
            print('   ... and %d more' % (len(degenerate) - 8))
    if not rows:
        raise SystemExit('\nnothing scorable: every combination was degenerate or had < %d shared drugs. '
                         'Check that the two cells\' label matrices actually contain the set members.'
                         % args.min_drugs)
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
    if moa_re is not None:
        ok = [r for r in rows if np.isfinite(r['r_moa'])]
        if ok:
            ok.sort(key=lambda r: r['p_moa'])
            qm = bh([r['p_moa'] for r in ok])
            for r, qq in zip(ok, qm):
                r['bh_q_moa'] = qq
            print()
            print('MoA-STRATIFIED: matching drugs (where the hypothesis applies) vs the complement (control)')
            h = ('%-9s %-20s %-11s %5s %8s %8s %6s %8s'
                 % ('pair', 'set', 'variant', 'n_moa', 'r_moa', 'BHq_moa', 'n_oth', 'r_oth'))
            print(h)
            print('-' * len(h))
            for r in ok[:args.top]:
                print('%-9s %-20s %-11s %5d %+8.3f %8.3g %6d %+8.3f'
                      % ('%s/%s' % (r['cell_x'], r['cell_y']), r['set'][:20], r['variant'],
                         r['n_moa'], r['r_moa'], r['bh_q_moa'], r['n_other'], r['r_other']))
            pm = sum(1 for r in ok if r['r_moa'] > 0)
            po = sum(1 for r in ok if np.isfinite(r['r_other']) and r['r_other'] > 0)
            print()
            print('predicted POSITIVE sign: %d of %d (%.0f%%) in the MATCHING group, '
                  '%d of %d (%.0f%%) in the COMPLEMENT'
                  % (pm, len(ok), 100 * pm / len(ok), po, len(ok), 100 * po / len(ok)))
            print('BH q < 0.05 in the matching group: %d, of which %d positive'
                  % (sum(1 for r in ok if r['bh_q_moa'] < 0.05),
                     sum(1 for r in ok if r['bh_q_moa'] < 0.05 and r['r_moa'] > 0)))
            print('median r: matching %+.3f, complement %+.3f  '
                  '(the hypothesis predicts the first to exceed the second)'
                  % (np.nanmedian([r['r_moa'] for r in ok]),
                     np.nanmedian([r['r_other'] for r in ok])))
    pos = sum(1 for r in rows if r['r_diff'] > 0)
    sig = [r for r in rows if r['bh_q'] < 0.05]
    print()
    print('%d of %d combinations have the predicted POSITIVE sign (%.0f%%; 50%% = no signal)'
          % (pos, len(rows), 100 * pos / len(rows)))
    print('%d pass BH q < 0.05, of which %d have the predicted positive sign'
          % (len(sig), sum(1 for r in sig if r['r_diff'] > 0)))
    print('median |naive_r| %.3f vs median |r_diff| %.3f - the naive value is the confounded one'
          % (np.nanmedian([abs(r['naive_r_in_x']) for r in rows]),
             np.nanmedian([abs(r['r_diff']) for r in rows])))
    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)) or '.', exist_ok=True)
        with open(args.out, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
            w.writeheader()
            w.writerows(rows)
        print('\nwrote %s' % args.out, file=sys.stderr)


if __name__ == '__main__':
    main()
