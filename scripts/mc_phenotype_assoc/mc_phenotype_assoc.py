#!/usr/bin/env python3
"""Do the MetChange labels separate drugs by a measured mitochondrial phenotype?

This is a LABEL-level test (link b of report 05, 9.3), not a model test: it asks
whether a drug's MetChange signature, summarised over a pathway's metabolites,
differs between drugs measured as active and inactive in an external assay
(Attene-Ramos membrane potential, Hallinger OCR).

Per drug statistic, two of them, both computed per perturbation and then averaged
over the drug's perturbations (doses of one drug are not independent):

  z      the set's mean log2FC standardised by that perturbation's own mean and
         SD over all features - how far the set sits from the drug's typical
         feature, in that drug's own units;
  auroc  the probability that a random set member outranks a random background
         feature in that perturbation - rank-based, so it is unmoved by the heavy
         tails of limma log2FC on rank-transformed MetChange scores.

The background for `auroc` is every feature that is not in the set under test.
Both statistics reduce to a mean over the set's rows of a per-row, per-column
quantity - the value itself for `z`, and the background ECDF of the value for
`auroc` - so the ECDF is computed once for the whole matrix and every set, real
or random, is then scored by one `nanmean`. That is what makes the random-set
control affordable.

Then, per set, the two drug groups are compared: Mann-Whitney U (two-sided),
group AUROC with a bootstrap CI, a label permutation p, and - the control that
matters most - an empirical p against `--n-random-sets` size-matched random
metabolite sets put through the identical pipeline. A pathway set has to beat
random sets of its own size, not just chance.

NO EXPECTED SIGN IS IMPOSED. A drug blocking the ETC post-translationally may
upregulate OXPHOS transcripts, so MetChange can read "more production" while flux
falls (the compensation hypothesis, report 05, 9.2). The test is two-sided and the
direction found is reported.

The activity filter is deliberately NOT the model's. Only ~14% of HepG2
perturbations have a significant MetChange feature, and requiring one costs most
of the power (147 vs 480 drugs becomes 45 vs 44). Since this test reads
continuous log2FC rather than a selected feature list, it runs on everything by
default; `--min-nsig 1` reruns it on the responsive subset. Run BOTH and compare:
if the effect is real it should be at least as strong in the responsive subset,
and that consistency is more convincing than either number alone.

Outputs in --out-dir:
  per_drug_scores.tsv   drug x set: z, auroc, n_perturbations, max_nsig
  set_results.tsv       set x contrast: group medians, AUROC + CI, p values, BH q
  set_member_values.npz compact (set members + background) x perturbation matrix,
                        so GSEA/NES or any other statistic can be added later
                        without re-reading the 1.5 GB label matrix
  run_info.json
"""
import argparse, collections, csv, json, os, sys, time

import numpy as np

try:
    from scipy import stats as sps
except Exception:                                      # noqa: BLE001
    sps = None


# ----------------------------------------------------------------- statistics
def mannwhitney(a, b):
    """Two-sided Mann-Whitney U with a tie-corrected normal approximation."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    a, b = a[np.isfinite(a)], b[np.isfinite(b)]
    n1, n2 = len(a), len(b)
    if n1 < 3 or n2 < 3:
        return np.nan, np.nan, np.nan
    if sps is not None:
        u, p = sps.mannwhitneyu(a, b, alternative='two-sided')
        return u, p, u / (n1 * n2)
    allv = np.concatenate([a, b])
    order = allv.argsort(kind='mergesort')
    ranks = np.empty(len(allv), float)
    ranks[order] = np.arange(1, len(allv) + 1)
    # average ranks within ties
    srt = allv[order]
    i = 0
    while i < len(srt):
        j = i
        while j + 1 < len(srt) and srt[j + 1] == srt[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = ranks[order[i:j + 1]].mean()
        i = j + 1
    r1 = ranks[:n1].sum()
    u = r1 - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    _, cnt = np.unique(allv, return_counts=True)
    n = n1 + n2
    tie = (cnt ** 3 - cnt).sum()
    var = n1 * n2 / 12.0 * ((n + 1) - tie / float(n * (n - 1)))
    if var <= 0:
        return u, np.nan, u / (n1 * n2)
    z = (abs(u - mu) - 0.5) / np.sqrt(var)
    from math import erfc, sqrt
    return u, erfc(z / sqrt(2.0)), u / (n1 * n2)


def auroc(case, ctrl):
    _, _, a = mannwhitney(case, ctrl)
    return a


def boot_ci(case, ctrl, n=2000, seed=0):
    rng = np.random.default_rng(seed)
    case, ctrl = np.asarray(case, float), np.asarray(ctrl, float)
    case, ctrl = case[np.isfinite(case)], ctrl[np.isfinite(ctrl)]
    if len(case) < 3 or len(ctrl) < 3:
        return np.nan, np.nan
    vals = [auroc(rng.choice(case, len(case), replace=True),
                  rng.choice(ctrl, len(ctrl), replace=True)) for _ in range(n)]
    return float(np.nanpercentile(vals, 2.5)), float(np.nanpercentile(vals, 97.5))


def perm_p(case, ctrl, n=5000, seed=0):
    """Permutation p for |AUROC - 0.5| with the group labels shuffled."""
    rng = np.random.default_rng(seed)
    case, ctrl = np.asarray(case, float), np.asarray(ctrl, float)
    case, ctrl = case[np.isfinite(case)], ctrl[np.isfinite(ctrl)]
    if len(case) < 3 or len(ctrl) < 3:
        return np.nan
    obs = abs(auroc(case, ctrl) - 0.5)
    pool = np.concatenate([case, ctrl])
    k = len(case)
    hits = 0
    for _ in range(n):
        rng.shuffle(pool)
        if abs(auroc(pool[:k], pool[k:]) - 0.5) >= obs:
            hits += 1
    return (hits + 1) / float(n + 1)


def bh(ps):
    ps = np.asarray(ps, float)
    ok = np.isfinite(ps)
    q = np.full(len(ps), np.nan)
    p = ps[ok]
    if not len(p):
        return q
    order = p.argsort()
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(len(p))
    out[order] = np.clip(ranked, 0, 1)
    q[ok] = out
    return q


# ----------------------------------------------------------------- data input
def read_sets(path, min_size):
    sets = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            members = [m for m in r['members'].split(',') if m]
            if len(members) >= min_size:
                sets.append((r['set'], r['variant'], set(members)))
    return sets


def read_phenotype(path):
    rows = list(csv.DictReader(open(path), delimiter='\t'))
    return {r['drug']: r for r in rows}


def read_fit_info(path):
    """perturbation -> (drug, n_sig)."""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            try:
                ns = int(r.get('n_sig') or 0)
            except ValueError:
                ns = 0
            out[r['perturbation']] = (r.get('drug', ''), ns)
    return out


def stream_matrix(path, keep_cols, report_every=500):
    """One pass over the feature x perturbation matrix, selected columns only.

    Keeps every feature row: 5,707 features x a few thousand perturbations is
    ~60 MB as float32, so there is no reason to subsample the background.
    """
    with open(path) as fh:
        header = next(fh).rstrip('\n').split(',')
        cols = header[1:]
        idx = [i for i, c in enumerate(cols) if c in keep_cols]
        sel = np.array(idx, dtype=np.int64)
        names = [cols[i] for i in idx]
        n = len(sel)
        print('matrix: %d columns, %d selected' % (len(cols), n), file=sys.stderr)
        if not n:
            raise SystemExit('no label columns match the phenotype drugs')

        feats, rows = [], []
        t0 = 0.0
        start = time.time()
        for line in fh:
            cut = line.find(',')
            feats.append(line[:cut].strip().strip('"'))
            parts = line[cut + 1:].rstrip('\n').split(',')
            vals = np.full(n, np.nan, dtype=np.float32)
            for k, i in enumerate(sel):
                try:
                    vals[k] = float(parts[i])
                except (IndexError, ValueError):
                    pass
            rows.append(vals)
            if len(feats) % report_every == 0:
                t0 = time.time() - start
                print('  %d rows (%.0fs)' % (len(feats), t0), file=sys.stderr)
    values = np.asarray(rows)
    print('  %d rows total (%.0fs)' % (len(feats), time.time() - start), file=sys.stderr)
    return names, feats, values


def column_stats(values):
    """Per-column mean / SD / n over all features, and the background ECDF of
    every cell: ecdf[i, j] = fraction of column j's other features below
    values[i, j], mid-rank for ties. AUROC of a set against the rest is then the
    mean of ecdf over the set's rows (exact up to the set's own members being
    included in the background, which is a 1/n_features effect)."""
    with np.errstate(invalid='ignore'):
        mean = np.nanmean(values, axis=0)
        sd = np.nanstd(values, axis=0)
    cnt = np.sum(np.isfinite(values), axis=0)
    ecdf = np.full(values.shape, np.nan, dtype=np.float32)
    for j in range(values.shape[1]):
        col = values[:, j]
        fin = np.isfinite(col)
        v = col[fin]
        if v.size < 2:
            continue
        order = v.argsort(kind='mergesort')
        ranks = np.empty(v.size, np.float64)
        ranks[order] = np.arange(1, v.size + 1)
        srt = v[order]
        i = 0
        while i < srt.size:
            k = i
            while k + 1 < srt.size and srt[k + 1] == srt[i]:
                k += 1
            if k > i:
                ranks[order[i:k + 1]] = ranks[order[i:k + 1]].mean()
            i = k + 1
        ecdf[fin, j] = ((ranks - 0.5) / v.size).astype(np.float32)
    return mean, sd, cnt, ecdf


# ----------------------------------------------------------------- the test
def score_sets(sets, feat_index, values, ecdf, mean, sd):
    """set x column -> (z, auroc), both as a mean over the set's rows."""
    out = {}
    for name, variant, members in sets:
        rows = [feat_index[m] for m in members if m in feat_index]
        if len(rows) < 3:
            continue
        with np.errstate(invalid='ignore'):
            m = np.nanmean(values[rows], axis=0)
            z = (m - mean) / np.where(sd > 0, sd, np.nan)
            a = np.nanmean(ecdf[rows], axis=0)
        out[(name, variant)] = (z, a)
    return out


def main():
    repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--labels', required=True, help='label dir with log2FC.csv + fit_info_nsig.csv')
    ap.add_argument('--matrix', default='log2FC.csv')
    ap.add_argument('--sets', default=os.path.join(repo, 'data/validation/metabolite_sets/metabolite_sets.tsv'))
    ap.add_argument('--phenotype', required=True)
    ap.add_argument('--case', default='MMP_down')
    ap.add_argument('--control', default='inactive')
    ap.add_argument('--extra-contrast', action='append', default=[],
                    help='further CASE:CONTROL pairs, e.g. MMP_up:inactive')
    ap.add_argument('--min-nsig', type=int, default=0,
                    help='0 = use every perturbation (default); 1 = MetChange-responsive only')
    ap.add_argument('--min-set-size', type=int, default=5)
    ap.add_argument('--n-random-sets', type=int, default=200)
    ap.add_argument('--n-perm', type=int, default=5000)
    ap.add_argument('--n-boot', type=int, default=2000)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    pheno = read_phenotype(args.phenotype)
    sets = read_sets(args.sets, args.min_set_size)
    fit = read_fit_info(os.path.join(args.labels, 'fit_info_nsig.csv'))
    contrasts = [(args.case, args.control)] + [tuple(c.split(':', 1)) for c in args.extra_contrast]
    wanted_groups = {g for c in contrasts for g in c}

    keep_cols, col_drug, col_nsig = set(), {}, {}
    for pert, (drug, nsig) in fit.items():
        p = pheno.get(drug)
        if not p or p['group'] not in wanted_groups:
            continue
        if nsig < args.min_nsig:
            continue
        keep_cols.add(pert)
        col_drug[pert] = drug
        col_nsig[pert] = nsig
    print('perturbations selected: %d (%d drugs)'
          % (len(keep_cols), len({col_drug[c] for c in keep_cols})), file=sys.stderr)

    need_rows = set()
    for _, _, m in sets:
        need_rows |= m

    names, feats, values = stream_matrix(
        os.path.join(args.labels, args.matrix), keep_cols)
    mean, sd, cnt, ecdf = column_stats(values)
    feat_index = {f: i for i, f in enumerate(feats)}
    print('features: %d (%d of them are in some set)'
          % (len(feats), sum(1 for f in feats if f in need_rows)), file=sys.stderr)

    scored = score_sets(sets, feat_index, values, ecdf, mean, sd)
    set_size = {(n, v): len([m for m in s_ if m in feat_index]) for n, v, s_ in sets}

    # random size-matched sets, identical pipeline: a pathway set must beat
    # random sets of its own size, not merely chance
    pool = [f for f in feats if f not in need_rows]
    sizes = sorted({set_size[k] for k in scored})
    rand_sets = []
    for size in sizes:
        if len(pool) <= size:
            continue
        for k in range(args.n_random_sets):
            rand_sets.append(('__random%d' % size, 'r%d' % k,
                              set(rng.choice(pool, size, replace=False))))
    rand_scored = score_sets(rand_sets, feat_index, values, ecdf, mean, sd) if rand_sets else {}
    print('random control: %d sets over %d sizes' % (len(rand_scored), len(sizes)),
          file=sys.stderr)

    # ------------------------------------------------------- aggregate per drug
    col_of = {c: i for i, c in enumerate(names)}
    drugs = sorted({col_drug[c] for c in names if c in col_drug})
    cols_by_drug = collections.defaultdict(list)
    for c in names:
        if c in col_drug:
            cols_by_drug[col_drug[c]].append(col_of[c])

    def per_drug(vec):
        return np.array([np.nanmean(vec[cols_by_drug[d]]) if cols_by_drug[d] else np.nan
                         for d in drugs])

    rows_out = []
    drug_stat = {}
    for (sname, variant), (z, a) in scored.items():
        dz, da = per_drug(z), per_drug(a)
        drug_stat[(sname, variant)] = (dz, da)
        for i, d in enumerate(drugs):
            rows_out.append([d, pheno[d]['group'], pheno[d].get('subgroup', ''), sname, variant,
                             '%.5g' % dz[i] if np.isfinite(dz[i]) else '',
                             '%.5g' % da[i] if np.isfinite(da[i]) else '',
                             len(cols_by_drug[d]),
                             max((col_nsig[c] for c in names
                                  if c in col_drug and col_drug[c] == d), default=0)])
    with open(os.path.join(args.out_dir, 'per_drug_scores.tsv'), 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['drug', 'group', 'subgroup', 'set', 'variant', 'z', 'auroc',
                    'n_perturbations', 'max_n_sig'])
        w.writerows(rows_out)

    rand_drug = {}
    for (sname, variant), (z, a) in rand_scored.items():
        rand_drug[(sname, variant)] = (per_drug(z), per_drug(a))

    # ------------------------------------------------------------ group tests
    groups = np.array([pheno[d]['group'] for d in drugs])
    results = []
    for case, ctrl in contrasts:
        ci_, cc_ = groups == case, groups == ctrl
        for (sname, variant), (dz, da) in sorted(drug_stat.items()):
            for stat_name, vec in (('z', dz), ('auroc', da)):
                a_, b_ = vec[ci_], vec[cc_]
                a_, b_ = a_[np.isfinite(a_)], b_[np.isfinite(b_)]
                if len(a_) < 3 or len(b_) < 3:
                    continue
                u, p, A = mannwhitney(a_, b_)
                lo, hi = boot_ci(a_, b_, args.n_boot, args.seed)
                pp = perm_p(a_, b_, args.n_perm, args.seed)
                # random size-matched sets: how often do they separate as well?
                size = set_size[(sname, variant)]
                rp, n_rand = np.nan, 0
                cand = [(k, v) for k, v in rand_drug.items() if k[0] == '__random%d' % size]
                if cand:
                    obs = abs(A - 0.5)
                    worse = 0
                    for _, (rz, ra) in cand:
                        rv = rz if stat_name == 'z' else ra
                        x_, y_ = rv[ci_], rv[cc_]
                        x_, y_ = x_[np.isfinite(x_)], y_[np.isfinite(y_)]
                        if len(x_) < 3 or len(y_) < 3:
                            continue
                        n_rand += 1
                        if abs(auroc(x_, y_) - 0.5) >= obs:
                            worse += 1
                    if n_rand:
                        rp = (worse + 1) / float(n_rand + 1)
                results.append({
                    'contrast': '%s_vs_%s' % (case, ctrl), 'set': sname, 'variant': variant,
                    'statistic': stat_name, 'n_members': size,
                    'n_case': len(a_), 'n_control': len(b_),
                    'median_case': float(np.median(a_)), 'median_control': float(np.median(b_)),
                    'delta_median': float(np.median(a_) - np.median(b_)),
                    'group_auroc': float(A), 'auroc_lo': lo, 'auroc_hi': hi,
                    'mannwhitney_p': float(p), 'perm_p': pp, 'random_set_p': rp,
                    'n_random_sets': n_rand,
                })
    for case, ctrl in contrasts:
        for stat_name in ('z', 'auroc'):
            sel = [r for r in results
                   if r['contrast'] == '%s_vs_%s' % (case, ctrl) and r['statistic'] == stat_name]
            qs = bh([r['mannwhitney_p'] for r in sel])
            for r, q in zip(sel, qs):
                r['bh_q'] = float(q)
    results.sort(key=lambda r: (r['contrast'], r['statistic'], r['mannwhitney_p']))
    cols = ['contrast', 'set', 'variant', 'statistic', 'n_members', 'n_case', 'n_control',
            'median_case', 'median_control', 'delta_median', 'group_auroc', 'auroc_lo',
            'auroc_hi', 'mannwhitney_p', 'perm_p', 'random_set_p', 'n_random_sets', 'bh_q']
    with open(os.path.join(args.out_dir, 'set_results.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(results)

    # keep the set members (small) plus the column annotation, so GSEA / NES or any
    # other statistic can be added later without re-reading the 1.5 GB matrix
    keep = sorted({feat_index[m] for _, _, mm in sets for m in mm if m in feat_index})
    np.savez_compressed(os.path.join(args.out_dir, 'set_member_values.npz'),
                        features=np.array([feats[i] for i in keep]),
                        columns=np.array(names), values=values[keep],
                        ecdf=ecdf[keep], col_mean=mean, col_sd=sd, col_n=cnt,
                        col_drug=np.array([col_drug.get(c, '') for c in names]),
                        col_nsig=np.array([col_nsig.get(c, 0) for c in names]))
    with open(os.path.join(args.out_dir, 'run_info.json'), 'w') as fh:
        json.dump({'labels': args.labels, 'phenotype': args.phenotype, 'sets': args.sets,
                   'contrasts': ['%s_vs_%s' % c for c in contrasts],
                   'min_nsig': args.min_nsig, 'n_features': len(feats),
                   'n_perturbations': len(names), 'n_drugs': len(drugs),
                   'n_sets_scored': len(scored), 'n_random_sets_per_size': args.n_random_sets,
                   'seed': args.seed,
                   'scipy': sps is not None}, fh, indent=2)

    print('\ntop rows of set_results.tsv:', file=sys.stderr)
    for r in results[:12]:
        print('  %-18s %-10s %-6s n=%d/%d  AUROC %.3f [%.3f,%.3f]  p=%.2g perm=%.3g rand=%.3g'
              % (r['set'], r['variant'], r['statistic'], r['n_case'], r['n_control'],
                 r['group_auroc'], r['auroc_lo'], r['auroc_hi'], r['mannwhitney_p'],
                 r['perm_p'], r['random_set_p']), file=sys.stderr)


if __name__ == '__main__':
    main()
