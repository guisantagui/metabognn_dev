#!/usr/bin/env python3
"""Do the MetChange labels separate drugs by a measured mitochondrial phenotype?

This is a LABEL-level test (link b of report 05, 9.3), not a model test: it asks
whether a drug's MetChange signature, summarised over a pathway's metabolites,
differs between drugs measured as active and inactive in an external assay
(Attene-Ramos membrane potential, Hallinger OCR).

Three modalities, selected with `--mode` (repeatable). Every statistic is computed
per perturbation and then averaged over the drug's perturbations, because doses of
one drug are not independent.

  --mode value   (magnitude)
    z             the set's mean log2FC standardised by that perturbation's own
                  mean and SD over all features - how far the set sits from the
                  drug's typical feature, in that drug's own units;
    auroc         the probability that a random set member outranks a random
                  background feature in that perturbation - rank-based, so it is
                  unmoved by the heavy tails of limma log2FC.

  --mode sign    (direction only, magnitude discarded)
    sign_mean     mean of sign(log2FC) over the set, in [-1, 1]: +1 = every member
                  up, -1 = every member down;
    sign_dev      the same minus the perturbation's own mean sign, so a drug that
                  pushes everything one way does not score on every set. This is
                  the directional reading: it asks whether the pathway moves UP or
                  DOWN, where `value` asks whether it moves at all.

  --mode sig     (significance, needs --qmatrix; USE WITH --min-nsig 1)
    sig_dev       fraction of the set that is significant at --qvalue, minus the
                  perturbation's overall significant fraction: an enrichment of
                  responding features in the pathway, the per-drug analogue of
                  `scripts/schuhknecht_mc_assoc/overlap_enrichment.py`;
    sig_signed_dev  mean of sign(log2FC) x 1[q < threshold] over the set, minus the
                  column mean: direction AND significance together, so features
                  that merely wobble contribute nothing.

**The `sig` statistics are only meaningful on perturbations that have significant
features at all.** With no significant feature the set fraction and the column
fraction are both 0, so `sig_dev` is identically 0; pooling those drugs with
responsive ones makes the test detect "has any significant feature" rather than
"the pathway is enriched". Measured on the HepG2 expression labels: the
non-responsive stratum gives AUROC exactly 0.500 (p = 1, all ties), while
responsive vs non-responsive drugs *within the controls* separate at AUROC 0.087
(p = 6e-60). So run `--mode sig` together with `--min-nsig 1`, and read it as a
statement about responsive drugs only. `value` and `sign` do not have this
problem: they are defined for every perturbation.

`value` and `sign` answer different questions and can disagree: a set can shift in
magnitude with mixed signs (value yes, sign no), or move consistently but weakly
(sign yes, value no). Report whichever were run, not the best one.

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
  per_drug_scores.tsv   drug x set x statistic, with n_perturbations and max_nsig
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
    """Returns the sets with identical member lists collapsed, plus the aliases.

    Sets built from the same subsystem can come out identical: no gene of these
    subsystems passes the hub cutoff, so every `*_nohub` gene set equals its
    `gene` set, and on the metabolite side `glycolysis_endpoints all_nohub`
    equals `all`. Testing the same members twice adds no information but doubles
    the BH family and so inflates every q-value, which is why they are collapsed
    here rather than merely flagged. The dropped names are kept in `aliases` so
    nothing is silently lost.
    """
    raw = []
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            members = [m for m in r['members'].split(',') if m]
            if len(members) >= min_size:
                raw.append((r['set'], r['variant'], set(members)))
    by_members = collections.OrderedDict()
    for name, variant, members in raw:
        by_members.setdefault(frozenset(members), []).append((name, variant))
    sets, aliases = [], {}
    for members, labels in by_members.items():
        # keep the plainest variant name: no "_nohub" suffix if one exists
        labels.sort(key=lambda nv: ('_nohub' in nv[1], nv[1], nv[0]))
        name, variant = labels[0]
        sets.append((name, variant, set(members)))
        aliases[(name, variant)] = ';'.join('%s %s' % nv for nv in labels[1:])
    n_dup = len(raw) - len(sets)
    if n_dup:
        print('sets: %d rows -> %d distinct member lists (%d duplicate%s collapsed)'
              % (len(raw), len(sets), n_dup, '' if n_dup == 1 else 's'), file=sys.stderr)
        for k, v in aliases.items():
            if v:
                print('   %s %s  ==  %s' % (k[0], k[1], v), file=sys.stderr)
    return sets, aliases


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


def _parse_pert_name(pert):
    """'vorinostat_0.66 uM_24 h' -> ('vorinostat', '0.66', '24'); used when a label
    directory has no fit_info_nsig.csv (the expression dirs may not)."""
    bits = pert.rsplit('_', 2)
    if len(bits) != 3:
        return pert, '', ''
    drug, dose, time = bits
    return drug, dose.split()[0] if dose.split() else '', time.split()[0] if time.split() else ''


def stream_matrices(paths, keep_cols, report_every=1000):
    """One pass over one or two aligned feature x perturbation matrices.

    `paths` maps a name ('value', 'q') to a file. The files are written by the same
    label pipeline, so they must share the header and the row order; both are
    checked as we go rather than assumed, because silently misaligning q-values
    with fold changes would invent significance.
    """
    handles = {k: open(v) for k, v in paths.items()}
    try:
        headers = {k: next(fh).rstrip('\n').split(',') for k, fh in handles.items()}
        ref = headers[list(paths)[0]]
        for k, h in headers.items():
            if h != ref:
                raise SystemExit('matrix %s has a different header to %s'
                                 % (k, list(paths)[0]))
        cols = ref[1:]
        idx = [i for i, c in enumerate(cols) if c in keep_cols]
        sel = np.array(idx, dtype=np.int64)
        names = [cols[i] for i in idx]
        n = len(sel)
        print('matrix: %d columns, %d selected' % (len(cols), n), file=sys.stderr)
        if not n:
            raise SystemExit('no label columns match the phenotype drugs')

        feats = []
        out = {k: [] for k in paths}
        start_t = time.time()
        while True:
            lines = {}
            for k, fh in handles.items():
                lines[k] = fh.readline()
            if not lines[list(paths)[0]]:
                break
            row_ids = set()
            for k, line in lines.items():
                if not line:
                    raise SystemExit('matrix %s ended early' % k)
                cut = line.find(',')
                row_ids.add(line[:cut].strip().strip('"'))
                parts = line[cut + 1:].rstrip('\n').split(',')
                vals = np.full(n, np.nan, dtype=np.float32)
                for kk, i in enumerate(sel):
                    try:
                        vals[kk] = float(parts[i])
                    except (IndexError, ValueError):
                        pass
                out[k].append(vals)
            if len(row_ids) != 1:
                raise SystemExit('row order differs between matrices: %s' % row_ids)
            feats.append(row_ids.pop())
            if len(feats) % report_every == 0:
                print('  %d rows (%.0fs)' % (len(feats), time.time() - start_t),
                      file=sys.stderr)
    finally:
        for fh in handles.values():
            fh.close()
    print('  %d rows total (%.0fs)' % (len(feats), time.time() - start_t), file=sys.stderr)
    return names, feats, {k: np.asarray(v) for k, v in out.items()}


def within_column_ecdf(values):
    """ecdf[i, j] = (midrank of values[i, j] within column j) / n, so that the mean
    of ecdf over a set's rows IS that set's AUROC against the other features."""
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
    return ecdf


def build_cell_matrices(modes, mats, qvalue):
    """name -> (per-cell matrix, how to baseline it against the column).

    Every statistic in this script is a mean over the set's rows of one such
    matrix, optionally re-expressed against the whole column. Keeping them in one
    place is what lets `value`, `sign` and `sig` - and the size-matched random-set
    control - share a single code path.
      'z'   subtract the column mean, divide by the column SD
      'dev' subtract the column mean
      'raw' leave as is
    """
    cells = {}
    vals = mats.get('value')
    qs = mats.get('q')
    if 'value' in modes:
        if vals is None:
            raise SystemExit('--mode value needs --matrix')
        cells['z'] = (vals, 'z')
        cells['auroc'] = (within_column_ecdf(vals), 'raw')
    if 'sign' in modes:
        if vals is None:
            raise SystemExit('--mode sign needs --matrix')
        sg = np.sign(vals).astype(np.float32)
        sg[~np.isfinite(vals)] = np.nan
        cells['sign_mean'] = (sg, 'raw')
        cells['sign_dev'] = (sg, 'dev')
    if 'sig' in modes:
        if qs is None:
            raise SystemExit('--mode sig needs --qmatrix')
        sig = (qs < qvalue).astype(np.float32)
        sig[~np.isfinite(qs)] = np.nan
        cells['sig_dev'] = (sig, 'dev')
        if vals is not None:
            ss = sig * np.sign(vals).astype(np.float32)
            ss[~np.isfinite(vals) | ~np.isfinite(qs)] = np.nan
            cells['sig_signed_dev'] = (ss, 'dev')
    return cells


def column_baselines(cells):
    out = {}
    for name, (m, how) in cells.items():
        with np.errstate(invalid='ignore'):
            mean = np.nanmean(m, axis=0)
            sd = np.nanstd(m, axis=0)
        out[name] = (mean, sd, how)
    return out


def score_sets(sets, feat_index, cells, baselines):
    """(set, variant) -> {statistic: value per column}."""
    out = {}
    for name, variant, members in sets:
        rows = [feat_index[m] for m in members if m in feat_index]
        if len(rows) < 3:
            continue
        stats = {}
        for stat, (m, _) in cells.items():
            mean, sd, how = baselines[stat]
            with np.errstate(invalid='ignore'):
                v = np.nanmean(m[rows], axis=0)
                if how == 'z':
                    v = (v - mean) / np.where(sd > 0, sd, np.nan)
                elif how == 'dev':
                    v = v - mean
            stats[stat] = v
        out[(name, variant)] = stats
    return out


def main():
    repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--labels', required=True, help='label dir with log2FC.csv + fit_info_nsig.csv')
    ap.add_argument('--matrix', default='log2FC.csv',
                    help='fold-change matrix inside --labels; needed by modes value/sign')
    ap.add_argument('--qmatrix', default='qvals.csv',
                    help='q-value matrix inside --labels; needed by mode sig')
    ap.add_argument('--mode', action='append', default=[],
                    choices=['value', 'sign', 'sig'],
                    help='repeatable; default value')
    ap.add_argument('--qvalue', type=float, default=0.05)
    ap.add_argument('--no-nsig', action='store_true',
                    help='do not read the q-value matrix just to record n_sig; the '
                         'responsiveness stratification is then unavailable')
    ap.add_argument('--fit-info', default='fit_info_nsig.csv',
                    help='optional; if absent, drug/dose/time are parsed from the '
                         'column names and n_sig is counted from the q-value matrix')
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
    ap.add_argument('--random-set-sizes', default='',
                    help='extra absolute sizes to include in the random null, beyond the '
                         'real sets own sizes. Used to compare the SHAPE of the null '
                         'across feature spaces at matched universe fraction (e.g. 42 of '
                         '5,707 metabolites = 0.74%% matches 149 of 20,186 genes).')
    ap.add_argument('--dump-random-null', action='store_true',
                    help='write random_null.tsv: one row per random set with its group '
                         'AUROC, so the null distributions can be compared directly')
    ap.add_argument('--n-perm', type=int, default=5000)
    ap.add_argument('--n-boot', type=int, default=2000)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    rng = np.random.default_rng(args.seed)

    modes = args.mode or ['value']
    pheno = read_phenotype(args.phenotype)
    sets, set_aliases = read_sets(args.sets, args.min_set_size)
    fit_path = os.path.join(args.labels, args.fit_info)
    fit = read_fit_info(fit_path) if os.path.exists(fit_path) else None
    if fit is None:
        print('no %s: drug/dose/time from column names, n_sig from the q-values'
              % args.fit_info, file=sys.stderr)
    contrasts = [(args.case, args.control)] + [tuple(c.split(':', 1)) for c in args.extra_contrast]
    wanted_groups = {g for c in contrasts for g in c}

    def pert_drug(pert):
        return fit[pert][0] if fit is not None else _parse_pert_name(pert)[0]

    if fit is not None:
        all_perts = list(fit)
    else:
        with open(os.path.join(args.labels,
                               args.matrix if 'value' in modes else args.qmatrix)) as fh:
            all_perts = next(fh).rstrip('\n').split(',')[1:]

    keep_cols, col_drug = set(), {}
    for pert in all_perts:
        drug = pert_drug(pert)
        p = pheno.get(drug)
        if not p or p['group'] not in wanted_groups:
            continue
        if fit is not None and fit[pert][1] < args.min_nsig:
            continue
        keep_cols.add(pert)
        col_drug[pert] = drug
    print('perturbations selected: %d (%d drugs)'
          % (len(keep_cols), len({col_drug[c] for c in keep_cols})), file=sys.stderr)

    need_rows = set()
    for _, _, m in sets:
        need_rows |= m

    paths = {}
    if 'value' in modes or 'sign' in modes:
        paths['value'] = os.path.join(args.labels, args.matrix)
    qpath = os.path.join(args.labels, args.qmatrix)
    # The q-values are read when mode sig needs them, and ALSO whenever there is no
    # fit_info_nsig.csv to take n_sig from: without n_sig the responsiveness
    # stratification of stratified_check.py cannot run, and that check is how the
    # main confounder of this analysis is controlled. --no-nsig skips the extra pass.
    if 'sig' in modes or (fit is None and not args.no_nsig and os.path.exists(qpath)):
        paths['q'] = qpath
    if 'sig' not in modes and 'q' in paths:
        print('reading %s as well, for n_sig (use --no-nsig to skip)' % args.qmatrix,
              file=sys.stderr)
    names, feats, mats = stream_matrices(paths, keep_cols)

    # n_sig per column: from fit_info when available, else counted from the q-values
    if fit is not None:
        col_nsig = {c: fit[c][1] for c in names}
    elif 'q' in mats:
        counts = np.nansum(mats['q'] < args.qvalue, axis=0)
        col_nsig = {c: int(counts[i]) for i, c in enumerate(names)}
    else:
        col_nsig = {c: 0 for c in names}
    if fit is None and args.min_nsig > 0:
        drop = [i for i, c in enumerate(names) if col_nsig[c] < args.min_nsig]
        if drop:
            keep = np.array([i for i in range(len(names)) if i not in set(drop)])
            names = [names[i] for i in keep]
            mats = {k: v[:, keep] for k, v in mats.items()}
            col_nsig = {c: col_nsig[c] for c in names}
            print('min_nsig applied after reading: %d columns kept' % len(names),
                  file=sys.stderr)

    cells = build_cell_matrices(modes, mats, args.qvalue)
    baselines = column_baselines(cells)
    feat_index = {f: i for i, f in enumerate(feats)}
    need_rows = set()
    for _, _, m in sets:
        need_rows |= m
    print('features: %d (%d of them are in some set); statistics: %s'
          % (len(feats), sum(1 for f in feats if f in need_rows),
             ', '.join(sorted(cells))), file=sys.stderr)

    scored = score_sets(sets, feat_index, cells, baselines)
    set_size = {(n, v): len([m for m in s_ if m in feat_index]) for n, v, s_ in sets}

    # The sets are built from Human-GEM alone and are not restricted to any one
    # label matrix, so report what was actually found here: a set whose members
    # are largely absent from this cell line's matrix is scored on a fraction of
    # itself, and that has to be visible rather than silent.
    short = [(n, v, set_size[(n, v)], len(s_)) for n, v, s_ in sets
             if (n, v) in set_size and set_size[(n, v)] < len(s_)]
    if short:
        n_bad = sum(1 for _, _, f, _ in short if f < args.min_set_size)
        print('sets with members absent from this label matrix: %d of %d'
              % (len(short), len(sets)), file=sys.stderr)
        for n, v, f, tot in sorted(short, key=lambda r: r[2] / float(r[3]))[:10]:
            print('   %-28s %-11s %d of %d members found (%.0f%%)'
                  % (n, v, f, tot, 100.0 * f / tot), file=sys.stderr)
        if n_bad:
            print('   %d set(s) fall below --min-set-size here and are not scored'
                  % n_bad, file=sys.stderr)

    # random size-matched sets, identical pipeline: a pathway set must beat random
    # sets of its own size, not merely chance
    pool = [f for f in feats if f not in need_rows]
    sizes = sorted({set_size[k] for k in scored}
                   | {int(x) for x in args.random_set_sizes.split() if x.strip()})
    rand_sets = []
    for size in sizes:
        if len(pool) <= size:
            continue
        for k in range(args.n_random_sets):
            rand_sets.append(('__random%d' % size, 'r%d' % k,
                              set(rng.choice(pool, size, replace=False))))
    rand_scored = score_sets(rand_sets, feat_index, cells, baselines) if rand_sets else {}
    print('random control: %d sets over %d sizes' % (len(rand_scored), len(sizes)),
          file=sys.stderr)

    # ------------------------------------------------------- aggregate per drug
    col_of = {c: i for i, c in enumerate(names)}
    drugs = sorted({col_drug[c] for c in names if c in col_drug})
    cols_by_drug = collections.defaultdict(list)
    for c in names:
        if c in col_drug:
            cols_by_drug[col_drug[c]].append(col_of[c])

    drug_nsig = {d: max((col_nsig[c] for c in names
                         if col_drug.get(c) == d), default=0) for d in drugs}

    def per_drug(vec):
        return np.array([np.nanmean(vec[cols_by_drug[d]]) if cols_by_drug[d] else np.nan
                         for d in drugs])

    stat_names = sorted(cells)
    rows_out = []
    drug_stat = {}
    for (sname, variant), stats in scored.items():
        per = {st: per_drug(v) for st, v in stats.items()}
        drug_stat[(sname, variant)] = per
        for i, d in enumerate(drugs):
            rows_out.append([d, pheno[d]['group'], pheno[d].get('subgroup', ''), sname, variant,
                             len(cols_by_drug[d]), drug_nsig[d]]
                            + ['%.5g' % per[st][i] if np.isfinite(per[st][i]) else ''
                               for st in stat_names])
    with open(os.path.join(args.out_dir, 'per_drug_scores.tsv'), 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(['drug', 'group', 'subgroup', 'set', 'variant',
                    'n_perturbations', 'max_n_sig'] + stat_names)
        w.writerows(rows_out)

    rand_drug = {k: {st: per_drug(v) for st, v in stats.items()}
                 for k, stats in rand_scored.items()}

    if args.dump_random_null and rand_drug:
        with open(os.path.join(args.out_dir, 'random_null.tsv'), 'w', newline='') as fh:
            w = csv.writer(fh, delimiter='\t')
            w.writerow(['contrast', 'size', 'universe_fraction', 'replicate',
                        'statistic', 'group_auroc', 'abs_effect'])
            groups0 = np.array([pheno[d]['group'] for d in drugs])
            for case, ctrl in contrasts:
                ci0, cc0 = groups0 == case, groups0 == ctrl
                for (rname, rvar), per_r in sorted(rand_drug.items()):
                    size = int(rname.replace('__random', ''))
                    for st, vec in per_r.items():
                        x_, y_ = vec[ci0], vec[cc0]
                        x_, y_ = x_[np.isfinite(x_)], y_[np.isfinite(y_)]
                        if len(x_) < 3 or len(y_) < 3:
                            continue
                        a_ = auroc(x_, y_)
                        w.writerow(['%s_vs_%s' % (case, ctrl), size,
                                    '%.6g' % (size / float(len(feats))), rvar, st,
                                    '%.6g' % a_, '%.6g' % abs(a_ - 0.5)])
        print('wrote random_null.tsv', file=sys.stderr)

    # ------------------------------------------------------------ group tests
    groups = np.array([pheno[d]['group'] for d in drugs])
    results = []
    for case, ctrl in contrasts:
        ci_, cc_ = groups == case, groups == ctrl
        for (sname, variant), per in sorted(drug_stat.items()):
            for stat_name in stat_names:
                vec = per[stat_name]
                a_, b_ = vec[ci_], vec[cc_]
                a_, b_ = a_[np.isfinite(a_)], b_[np.isfinite(b_)]
                if len(a_) < 3 or len(b_) < 3:
                    continue
                u, p, A = mannwhitney(a_, b_)
                lo, hi = boot_ci(a_, b_, args.n_boot, args.seed)
                pp = perm_p(a_, b_, args.n_perm, args.seed)
                size = set_size[(sname, variant)]
                rp, n_rand = np.nan, 0
                cand = [v for k, v in rand_drug.items() if k[0] == '__random%d' % size]
                if cand:
                    obs = abs(A - 0.5)
                    worse = 0
                    for per_r in cand:
                        rv = per_r[stat_name]
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
                    'aliases': set_aliases.get((sname, variant), ''),
                    'statistic': stat_name, 'n_members': size,
                    'n_case': len(a_), 'n_control': len(b_),
                    'median_case': float(np.median(a_)), 'median_control': float(np.median(b_)),
                    'delta_median': float(np.median(a_) - np.median(b_)),
                    'group_auroc': float(A), 'auroc_lo': lo, 'auroc_hi': hi,
                    'mannwhitney_p': float(p), 'perm_p': pp, 'random_set_p': rp,
                    'n_random_sets': n_rand,
                })
    for case, ctrl in contrasts:
        for stat_name in stat_names:
            sel = [r for r in results
                   if r['contrast'] == '%s_vs_%s' % (case, ctrl) and r['statistic'] == stat_name]
            qs = bh([r['mannwhitney_p'] for r in sel])
            for r, q in zip(sel, qs):
                r['bh_q'] = float(q)
    results.sort(key=lambda r: (r['contrast'], r['statistic'], r['mannwhitney_p']))
    cols = ['contrast', 'set', 'variant', 'aliases', 'statistic', 'n_members', 'n_case', 'n_control',
            'median_case', 'median_control', 'delta_median', 'group_auroc', 'auroc_lo',
            'auroc_hi', 'mannwhitney_p', 'perm_p', 'random_set_p', 'n_random_sets', 'bh_q']
    with open(os.path.join(args.out_dir, 'set_results.tsv'), 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(results)

    # keep the set members (small) plus the column annotation, so GSEA / NES or any
    # other statistic can be added later without re-reading the 1.5 GB matrix
    keep = sorted({feat_index[m] for _, _, mm in sets for m in mm if m in feat_index})
    dump = {'features': np.array([feats[i] for i in keep]),
            'columns': np.array(names),
            'col_drug': np.array([col_drug.get(c, '') for c in names]),
            'col_nsig': np.array([col_nsig.get(c, 0) for c in names])}
    for st, (m, _) in cells.items():
        dump['cells_' + st] = m[keep]
        dump['colmean_' + st] = baselines[st][0]
        dump['colsd_' + st] = baselines[st][1]
    np.savez_compressed(os.path.join(args.out_dir, 'set_member_values.npz'), **dump)
    with open(os.path.join(args.out_dir, 'run_info.json'), 'w') as fh:
        json.dump({'labels': args.labels, 'phenotype': args.phenotype, 'sets': args.sets,
                   'modes': modes, 'statistics': stat_names, 'qvalue': args.qvalue,
                   'matrix': args.matrix if 'value' in paths else None,
                   'qmatrix': args.qmatrix if 'q' in paths else None,
                   'contrasts': ['%s_vs_%s' % c for c in contrasts],
                   'min_nsig': args.min_nsig, 'n_features': len(feats),
                   'n_perturbations': len(names), 'n_drugs': len(drugs),
                   'n_sets_scored': len(scored), 'n_sets_in_bh_family': len(scored),
                   'n_random_sets_per_size': args.n_random_sets,
                   'seed': args.seed,
                   'scipy': sps is not None}, fh, indent=2)

    print('\ntop rows of set_results.tsv:', file=sys.stderr)
    for r in sorted(results, key=lambda r: r['mannwhitney_p'])[:12]:
        print('  %-18s %-10s %-6s n=%d/%d  AUROC %.3f [%.3f,%.3f]  p=%.2g perm=%.3g rand=%.3g'
              % (r['set'], r['variant'], r['statistic'], r['n_case'], r['n_control'],
                 r['group_auroc'], r['auroc_lo'], r['auroc_hi'], r['mannwhitney_p'],
                 r['perm_p'], r['random_set_p']), file=sys.stderr)


if __name__ == '__main__':
    main()
