#!/usr/bin/env python3
"""Is the MetChange / phenotype association a chemistry effect or a metabolic one?

Two findings make this the decisive question rather than an optional extra:

  - the Attene-Ramos TCA signal is carried entirely by compounds NOT in the
    ToxCast library (AUROC 0.352, n = 109/425) and is absent in the ToxCast
    subset (0.447, n = 38/55) - which is also why Hallinger, which *is* the
    ToxCast library, is null;
  - matching on n_sig equalises how MUCH a drug responds, not WHAT KIND of
    compound it is.

Tests, from weakest assumption to strongest:

  moa-strat-fine    van Elteren stratified Wilcoxon within each MoA class, using
                    ONE class per drug (see --moa-pick) so no drug is counted
                    twice. Asks whether phenotype-active drugs score differently
                    from inactive drugs OF THE SAME MECHANISM.
  moa-strat-coarse  the same over ~9 coarse mechanism families, so the strata are
                    large enough to have power. The fine classes are 225 free-text
                    labels of wildly varying granularity ("Cyclooxygenase
                    inhibitor" against "RET tyrosine kinase inhibitor"), which
                    makes the fine test both low-powered and restricted to a
                    population of receptor-targeting drugs.
  sim-matched       each case paired with its most chemically similar control
                    (ECFP4 Tanimoto, greedy, no replacement) and tested PAIRED
                    (Wilcoxon signed-rank on the within-pair differences). The
                    achieved similarity is reported, because this control is
                    worthless if the pairs are not actually similar.
  potency           within actives only, is the per-drug score monotone in the
                    MEASURED potency of the phenotype? A gradient is far harder
                    for a chemotype to produce than a binary difference.

`--spec` may be repeated to run a panel in one call; p-values are then corrected
across the panel per test type (Holm and BH both reported), because running five
set/statistic combinations and quoting the best uncorrected p is not a result.

Reads per_drug_scores.tsv and the phenotype table. rdkit is needed only for
sim-matched.
"""
import argparse, collections, csv, os, re, sys

import numpy as np

try:
    from scipy import stats as sps
except Exception:                                      # noqa: BLE001
    sps = None

# Coarse mechanism families, checked in order: the first match wins, so the more
# specific patterns come first. Built from the MoA vocabulary actually present in
# compoundinfo_beta.txt, not from an external ontology.
COARSE = [
    ('kinase inhibitor', r'kinase inhibitor|tyrosine kinase|\bmTOR\b|\bPI3K\b|\bJAK\b|\bALK\b|'
                         r'\bEGFR\b|\bVEGFR\b|\bPDGFR\b|\bRAF\b|\bMEK\b|\bCDK\b|Aurora|\bGSK\b|'
                         r'\bSYK\b|\bBTK\b|\bAKT\b|\bFLT3\b|\bSRC\b|\bABL\b'),
    ('DNA / cytoskeleton / cytotoxic', r'topoisomerase|tubulin|\bDNA\b|HDAC|histone|proteasome|'
                                      r'ribonucleotide reductase|antimetabolite|alkylat|intercalat|'
                                      r'telomerase|PARP'),
    ('nuclear receptor ligand', r'glucocorticoid|estrogen|androgen|progest|retinoid|\bPPAR\b|'
                                r'thyroid|vitamin D|mineralocorticoid|\bRAR\b|\bRXR\b|'
                                r'nuclear receptor|steroid'),
    ('GPCR / neurotransmitter', r'dopamine|serotonin|adrenergic|histamine|acetylcholine|muscarinic|'
                                r'nicotinic|opioid|cannabinoid|\bGABA\b|glutamate|adenosine receptor|'
                                r'angiotensin|endothelin|prostaglandin receptor|leukotriene|'
                                r'bradykinin|vasopressin|oxytocin|melatonin'),
    ('ion channel', r'channel'),
    ('transporter / pump', r'transporter|reuptake|ATPase|\bpump\b|symporter|exchanger|'
                           r'proton pump|efflux'),
    ('protease / peptidase', r'protease|peptidase|convertase|thrombin|elastase'),
    ('metabolic enzyme', r'dehydrogenase|synthase|synthetase|reductase|oxidase|carboxylase|'
                         r'lipase|HMG|cyclooxygenase|lipoxygenase|phosphodiesterase|'
                         r'carbonic anhydrase|aromatase|hydroxylase|transferase|'
                         r'cytochrome P450|phosphatase|decarboxylase|isomerase|'
                         r'glycosylase|amidase|esterase'),
]


def coarse_moa(label):
    for name, pat in COARSE:
        if re.search(pat, label, re.I):
            return name
    return 'other'


# ------------------------------------------------------------------ statistics
def mw(a, b):
    """Unpaired Mann-Whitney; returns AUROC, p, n1, n2."""
    a = np.asarray([x for x in a if np.isfinite(x)], float)
    b = np.asarray([x for x in b if np.isfinite(x)], float)
    if len(a) < 1 or len(b) < 1 or sps is None:
        return np.nan, np.nan, len(a), len(b)
    u, p = sps.mannwhitneyu(a, b, alternative='two-sided')
    return u / (len(a) * len(b)), p, len(a), len(b)


def paired(case, ctrl):
    """Wilcoxon signed-rank on within-pair differences: the test the matched
    designs actually call for. Returns (fraction of pairs with case < control,
    p, n_pairs, median difference). The fraction is on the same 0-1 scale as an
    AUROC and means the same thing (below 0.5 = cases lower)."""
    d = np.array([c - k for c, k in zip(case, ctrl)], float)
    d = d[np.isfinite(d)]
    if len(d) < 6 or sps is None:
        return np.nan, np.nan, len(d), np.nan
    nz = d[d != 0]
    if len(nz) < 6:
        return np.nan, np.nan, len(d), float(np.median(d))
    _, p = sps.wilcoxon(nz, alternative='two-sided')
    return float(np.mean(d < 0)), p, len(d), float(np.median(d))


def van_elteren(strata):
    """Stratified Wilcoxon, weight 1/(n_s+1); mid-ranks, tie-corrected variance."""
    num = den = 0.0
    aurocs, weights, used = [], [], 0
    for a, b in strata:
        n1, n2 = len(a), len(b)
        if n1 < 1 or n2 < 1:
            continue
        n = n1 + n2
        allv = np.concatenate([a, b])
        order = allv.argsort(kind='mergesort')
        ranks = np.empty(n, float)
        ranks[order] = np.arange(1, n + 1)
        srt = allv[order]
        i = 0
        while i < n:
            k = i
            while k + 1 < n and srt[k + 1] == srt[i]:
                k += 1
            if k > i:
                ranks[order[i:k + 1]] = ranks[order[i:k + 1]].mean()
            i = k + 1
        w = 1.0 / (n + 1)
        r1 = ranks[:n1].sum()
        e = n1 * (n + 1) / 2.0
        _, cnt = np.unique(allv, return_counts=True)
        v = (n1 * n2 / 12.0 * ((n + 1) - (cnt ** 3 - cnt).sum() / float(n * (n - 1)))
             if n > 1 else 0.0)
        num += w * (r1 - e)
        den += (w ** 2) * v
        aurocs.append((r1 - n1 * (n1 + 1) / 2.0) / (n1 * n2))
        weights.append(w * n1 * n2)
        used += 1
    if den <= 0 or not used:
        return np.nan, np.nan, np.nan, used
    z = num / np.sqrt(den)
    from math import erfc, sqrt
    return z, erfc(abs(z) / sqrt(2.0)), float(np.average(aurocs, weights=weights)), used


def boot_strata(strata, n=2000, seed=0):
    """Cluster bootstrap over strata of the weighted AUROC."""
    rng = np.random.default_rng(seed)

    def wa(sts):
        num = den = 0.0
        for a, b in sts:
            n1, n2 = len(a), len(b)
            if not n1 or not n2:
                continue
            u = sum((x > y) + 0.5 * (x == y) for x in a for y in b)
            w = n1 * n2 / float(n1 + n2 + 1)
            num += w * (u / (n1 * n2))
            den += w
        return num / den if den else np.nan

    vals = [wa([strata[i] for i in rng.integers(0, len(strata), len(strata))])
            for _ in range(n)]
    return float(np.nanpercentile(vals, 2.5)), float(np.nanpercentile(vals, 97.5))


def holm_bh(ps):
    ps = np.asarray(ps, float)
    m = len(ps)
    o = ps.argsort()
    holm = np.empty(m)
    run = 0.0
    for rank, i in enumerate(o):
        run = max(run, (m - rank) * ps[i])
        holm[i] = min(run, 1.0)
    r = ps[o] * m / np.arange(1, m + 1)
    r = np.minimum.accumulate(r[::-1])[::-1]
    bh = np.empty(m)
    bh[o] = np.clip(r, 0, 1)
    return holm, bh


# ----------------------------------------------------------------------- input
def load_scores(path, set_name, variant, statistic):
    out = {}
    for r in csv.DictReader(open(path), delimiter='\t'):
        if r['set'] != set_name or r['variant'] != variant:
            continue
        try:
            v = float(r[statistic])
        except (ValueError, KeyError):
            v = np.nan
        out[r['drug']] = (r['group'], v)
    return out


def load_annotation(path, drugs):
    moa, smiles = collections.defaultdict(set), {}
    with open(path, encoding='utf-8', errors='replace') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            n = r.get('cmap_name', '')
            if n not in drugs:
                continue
            if r.get('moa'):
                moa[n].add(r['moa'].strip())
            s = (r.get('canonical_smiles') or '').strip()
            if s and s not in ('NA', 'restricted', '-666') and n not in smiles:
                smiles[n] = s
    return moa, smiles


def pick_primary(moa, how):
    """One MoA per drug, so no drug enters more than one stratum."""
    freq = collections.Counter(m for ms in moa.values() for m in ms)
    out = {}
    for d, ms in moa.items():
        if not ms:
            continue
        if how == 'common':
            out[d] = sorted(ms, key=lambda m: (-freq[m], m))[0]
        elif how == 'rare':
            out[d] = sorted(ms, key=lambda m: (freq[m], m))[0]
        else:
            out[d] = sorted(ms)[0]
    return out


def sim_match(cases, ctrls, smiles):
    from rdkit import Chem, DataStructs
    from rdkit.Chem import rdFingerprintGenerator
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    def fp(d):
        m = Chem.MolFromSmiles(smiles[d]) if d in smiles else None
        return gen.GetFingerprint(m) if m is not None else None

    cfp = {d: fp(d) for d in cases}
    kfp = {d: fp(d) for d in ctrls}
    cs = [d for d in cases if cfp[d] is not None]
    pool = [d for d in ctrls if kfp[d] is not None]
    order = sorted(cs, key=lambda d: -max(
        (DataStructs.TanimotoSimilarity(cfp[d], kfp[k]) for k in pool), default=0.0))
    pairs, sims = [], []
    for d in order:
        if not pool:
            break
        s = [DataStructs.TanimotoSimilarity(cfp[d], kfp[k]) for k in pool]
        j = int(np.argmax(s))
        pairs.append((d, pool[j]))
        sims.append(s[j])
        pool.pop(j)
    return pairs, np.array(sims)


# ------------------------------------------------------------------------ main
def run_one(spec, scores_path, pheno, moa, smiles, primary, args):
    set_name, variant, statistic, case = spec
    sc = load_scores(scores_path, set_name, variant, statistic)
    if not sc:
        print('  %s %s: set not found, skipped' % (set_name, variant), file=sys.stderr)
        return []
    val = {d: v for d, (_, v) in sc.items()}
    cases = [d for d, (g, v) in sc.items() if g == case and np.isfinite(v)]
    ctrls = [d for d, (g, v) in sc.items() if g == args.control and np.isfinite(v)]
    tag = '%s %s / %s / %s' % (set_name, variant, statistic, case)
    rows = []

    a, p, n1, n2 = mw([val[d] for d in cases], [val[d] for d in ctrls])
    rows.append(dict(spec=tag, test='reference (all drugs, unpaired)', effect=a, p=p,
                     n_case=n1, n_control=n2, note=''))

    for how, lab in (('fine', 'moa-strat-fine'), ('coarse', 'moa-strat-coarse')):
        key = (lambda d: primary.get(d)) if how == 'fine' else \
              (lambda d: coarse_moa(primary[d]) if d in primary else None)
        groups = collections.defaultdict(lambda: ([], []))
        for d in cases:
            k = key(d)
            if k:
                groups[k][0].append(val[d])
        for d in ctrls:
            k = key(d)
            if k:
                groups[k][1].append(val[d])
        strata = [(np.array(a_), np.array(b_)) for a_, b_ in groups.values() if a_ and b_]
        if not strata:
            continue
        z, pv, wa, used = van_elteren(strata)
        lo, hi = boot_strata(strata, args.n_boot, args.seed)
        nc = sum(len(s[0]) for s in strata)
        nk = sum(len(s[1]) for s in strata)
        rows.append(dict(spec=tag, test=lab, effect=wa, p=pv, n_case=nc, n_control=nk,
                         note='%d strata, CI [%.3f, %.3f], z=%.2f' % (used, lo, hi, z)))

    try:
        pairs, sims = sim_match(cases, ctrls, smiles)
        fr, pp, npair, md = paired([val[d] for d, _ in pairs], [val[k] for _, k in pairs])
        rows.append(dict(spec=tag, test='sim-matched (PAIRED signed-rank)', effect=fr, p=pp,
                         n_case=npair, n_control=npair,
                         note='median Tanimoto %.3f, median pair diff %.4g'
                              % (float(np.median(sims)), md)))
    except ImportError:
        pass

    pot, sco = [], []
    for d in cases:
        s = (pheno.get(d, {}).get('potency_uM') or '').strip()
        if s:
            try:
                pot.append(np.log10(float(s)))
                sco.append(val[d])
            except ValueError:
                pass
    if len(pot) > 10 and sps is not None:
        rho, prho = sps.spearmanr(pot, sco)
        rows.append(dict(spec=tag, test='potency gradient (Spearman, actives only)',
                         effect=rho, p=prho, n_case=len(pot), n_control=0,
                         note='rho vs log10 potency; + = more potent scores lower'))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--scores', required=True)
    ap.add_argument('--phenotype', required=True)
    ap.add_argument('--spec', action='append', required=True,
                    help='set:variant:statistic:case, repeatable')
    ap.add_argument('--control', default='inactive')
    ap.add_argument('--moa-pick', choices=['common', 'rare', 'alpha'], default='common',
                    help='which of a drugs MoA labels becomes its single stratum')
    ap.add_argument('--compoundinfo',
                    default='/Users/guillem.santamaria/Documents/postdoc/comput/'
                            'metabo_pert/data/LINCS/compoundinfo_beta.txt')
    ap.add_argument('--n-boot', type=int, default=2000)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--out')
    args = ap.parse_args()

    pheno = {r['drug']: r for r in csv.DictReader(open(args.phenotype), delimiter='\t')}
    all_drugs = {r['drug'] for r in csv.DictReader(open(args.scores), delimiter='\t')}
    moa, smiles = load_annotation(args.compoundinfo, all_drugs)
    primary = pick_primary(moa, args.moa_pick)
    multi = sum(1 for d, ms in moa.items() if len(ms) > 1)
    print('drugs in scores: %d | with MoA: %d (%d had >1 label, collapsed by --moa-pick %s) '
          '| fine classes: %d | coarse: %d'
          % (len(all_drugs), len(primary), multi, args.moa_pick,
             len(set(primary.values())), len({coarse_moa(m) for m in primary.values()})))
    print('coarse class sizes: %s' % dict(
        collections.Counter(coarse_moa(m) for m in primary.values()).most_common()))
    print()

    specs = [tuple(s.split(':')) for s in args.spec]
    rows = []
    for s in specs:
        if len(s) != 4:
            raise SystemExit('--spec must be set:variant:statistic:case, got %r' % (s,))
        rows += run_one(s, args.scores, pheno, moa, smiles, primary, args)

    # correction across the panel, per test type
    for test in sorted({r['test'] for r in rows}):
        sel = [r for r in rows if r['test'] == test]
        ps = [r['p'] for r in sel]
        if len(sel) > 1 and all(np.isfinite(ps)):
            holm, bh = holm_bh(ps)
            for r, h, b in zip(sel, holm, bh):
                r['p_holm'] = h
                r['p_bh'] = b
        else:
            for r in sel:
                r['p_holm'] = r['p']
                r['p_bh'] = r['p']

    w = '%-34s %-38s %8s %10s %10s %10s  %s'
    print(w % ('spec', 'test', 'effect', 'p', 'p_holm', 'p_BH', 'note'))
    print('-' * 150)
    for r in rows:
        print(w % (r['spec'][:34], r['test'][:38],
                   '%.3f' % r['effect'] if np.isfinite(r['effect']) else '-',
                   '%.3g' % r['p'] if np.isfinite(r['p']) else '-',
                   '%.3g' % r['p_holm'] if np.isfinite(r['p_holm']) else '-',
                   '%.3g' % r['p_bh'] if np.isfinite(r['p_bh']) else '-', r['note']))

    if args.out:
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        with open(args.out, 'w', newline='') as fh:
            ww = csv.DictWriter(fh, fieldnames=['spec', 'test', 'effect', 'p', 'p_holm',
                                                'p_bh', 'n_case', 'n_control', 'note'],
                                delimiter='\t')
            ww.writeheader()
            ww.writerows(rows)
        print('\nwrote %s' % args.out, file=sys.stderr)


if __name__ == '__main__':
    main()
