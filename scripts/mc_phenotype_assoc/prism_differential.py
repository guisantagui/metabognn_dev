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


def gem_gene_subsystems(nodes_path):
    """gene symbol -> {GEM subsystem}, from the reaction rows' GPRs."""
    out = {}
    with open(nodes_path, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r.get('node_type') != 'reaction' or not r.get('gpr') or not r.get('subsystem'):
                continue
            for g in re.split(r'[^A-Za-z0-9_.-]+', r['gpr']):
                if g and g not in ('or', 'and'):
                    out.setdefault(g, set()).add(r['subsystem'])
    return out


def drug_target_effects(path):
    """lowercased drug name -> {target symbol: effect}"""
    out = {}
    with open(path, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            n = (r.get('pertname') or '').strip().lower()
            if n:
                out.setdefault(n, {})[r['target']] = r.get('effect', '')
    return out


def jaccard_groups(vectors, threshold):
    """Single-linkage grouping of drugs by Jaccard similarity of their target sets.

    Union-find over pairs above the threshold. Drugs with no annotated targets each
    form their own group, so they are neither merged together nor discarded.
    """
    drugs = sorted(vectors)
    parent = {d: d for d in drugs}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra

    # only drugs sharing at least one target can exceed any positive threshold
    by_target = {}
    for d in drugs:
        for t in vectors[d]:
            by_target.setdefault(t, []).append(d)
    checked = set()
    for members in by_target.values():
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                a, b = members[i], members[j]
                key = (a, b) if a < b else (b, a)
                if key in checked:
                    continue
                checked.add(key)
                va, vb = vectors[a], vectors[b]
                if not va or not vb:
                    continue
                if len(va & vb) / len(va | vb) >= threshold:
                    union(a, b)
    groups = {}
    for d in drugs:
        groups.setdefault(find(d), []).append(d)
    return list(groups.values())


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
    ap.add_argument('--expr-run', default='',
                    help='name of the EXPRESSION run directory under --scores-root (normally '
                         '"prism_viability_expr"). When given, each row additionally reports the expression '
                         'differential and, more importantly, the PARTIAL correlation of the MetChange '
                         'differential with the phenotype CONTROLLING FOR the matching gene set. That is the '
                         'statistic the project turns on: the claim is not that MetChange beats the '
                         'transcriptome in a horse race - expression being informative too is expected - but '
                         'that integrating stoichiometry and the network carries information expression '
                         'alone does not. A MetChange association that vanishes once the gene set is '
                         'controlled for is a transcriptional association wearing a metabolic label.')
    ap.add_argument('--expr-set-map',
                    default='dntp_endpoints=nucleotide,ntp_endpoints=nucleotide,nmp_endpoints=nucleotide,'
                            'amino_acid_endpoints=amino_acid,tca_endpoints=tca,'
                            'glycolysis_endpoints=glycolysis,etc_endpoints=oxphos,'
                            'energy_currency=oxphos,adenylates=oxphos,redox_couples=oxphos',
                    help='which GENE set conditions which metabolite set, for the curated endpoint sets that '
                         'have no same-named gene counterpart (only the 10 subsystem sets do). Comma-separated '
                         'metabolite_set=gene_set pairs; anything not listed is conditioned on its own name.')
    ap.add_argument('--expr-variant', default='gene',
                    help="variant of the gene sets to condition on (default 'gene')")
    ap.add_argument('--targets', default='',
                    help='data/drug_targets/combined/drug_targets.tsv. Selects drugs by ANNOTATED TARGET '
                         'rather than by MoA free text, which is both principled and tells you which '
                         'subsystem should move. Needs --nodes for the gene -> subsystem map. Beware that '
                         '"has a GEM target" alone is far too permissive: 252 of 821 PRISM drugs qualify, but '
                         'the commonest subsystems are Transport (82 drugs), Arachidonic acid (69), Linoleate '
                         '(61), Leukotriene (58), Drug metabolism (54), Xenobiotics (52) and Inositol '
                         'phosphate (26, the PI3K inhibitors) - signalling lipids, xenobiotic handling and '
                         'transport, not biomass precursors. Hence --target-subsystems.')
    ap.add_argument('--nodes', default='',
                    help='data/networks/combined/nodes.tsv; reaction rows give gpr -> subsystem')
    ap.add_argument('--target-subsystems',
                    default=r'nucleotide|purine|pyrimidine|folate|one.carbon|amino acid|alanine|arginine|'
                            r'asparagine|aspartate|cysteine|glutam|glycine|histidine|isoleucine|leucine|'
                            r'lysine|methionine|phenylalanine|proline|serine|threonine|tryptophan|tyrosine|'
                            r'valine|urea cycle|glycolysis|citric acid|tricarboxylic|'
                            r'oxidative phosphorylation|pentose phosphate',
                    help='case-insensitive pattern on the GEM subsystem of a drug target. The default is the '
                         'biomass-precursor AND central-energy subsystems; with --target-effect inh that '
                         'selects 90 of the 821 PRISM drugs (49 if glycolysis/TCA/OXPHOS/PPP are dropped and '
                         'only the biomass precursors kept). Pass an empty string to accept any GEM target.')
    ap.add_argument('--relevance', default='',
                    help='target_relevance.tsv from build_target_relevance.py (gene, set, hops): the '
                         'shortest PRODUCING route from a gene to a member of each metabolite set. With it, a '
                         'drug is relevant to a set only if an inhibited target can actually affect that '
                         "set's production, which is derived from stoichiometry instead of from subsystem "
                         'labels. It excludes PDE, PI3K and EZH2 automatically (never reached) - the hand-made '
                         'list of report 10, 11.3 rederived. **Only use it with small endpoint sets**: for '
                         '`dntp_endpoints` it rules out 60% of GEM genes, but the subsystem sets hold 192-823 '
                         'metabolites so nearly every gene is within a few hops and it selects ~235 drugs at '
                         'any depth (report 10, 14).')
    ap.add_argument('--max-hops', type=int, default=4,
                    help='maximum producing-route distance for --relevance. Hop distance does NOT separate '
                         'biosynthetic routes from incidental ones - to a dNTP, TYMS and TK1 sit at 4 '
                         'alongside the G9a methyltransferase, DHFR and IMPDH1 at 5 alongside CYP2C19, DHODH '
                         'at 6 alongside MAO-A - so use --exclude-target-classes as well.')
    ap.add_argument('--min-target-sources', type=int, default=1,
                    help='minimum `n_sources` for a drug-target annotation to count. This is the fix for weak '
                         'annotations, which is a DIFFERENT failure from a misplaced enzyme: adaprev->HK1, '
                         'tofogliflozin->GCK and PJ-34->EEF2 all have n_sources=2 from CLUE/RepurposingHub, '
                         'while every genuine antimetabolite annotation carries 3-5 including ChEMBL and '
                         'DrugCentral. Those three genes were wrongly put in an exclusion list in report 10, '
                         '11.3; a hexokinase inhibitor would be a legitimate glycolysis drug. 3 is the '
                         'sensible value.')
    ap.add_argument('--exclude-target-classes', default='',
                    help='regex over target symbols for enzyme classes present in the GEM for catabolic, '
                         'cofactor-consuming, signalling or xenobiotic reasons rather than as biosynthetic '
                         'routes to a precursor. The documented default (report 10, 14.3) is '
                         "'^CYP|^UGT|^SULT|^NAT[12]$|^FMO|^EHMT|^EZH|^DNMT|^SETD|^KMT|^PRMT|^COMT$|^NNMT$|"
                         "^MAO|^ALDH|^ADH[0-9]|^AOX1$|^DBH$|^SIRT|^PARP|^CD38$|^TXNRD|^GSR$|^GPX|^GST|"
                         "^NQO1$|^SOD[0-9]|^CAT$|^PTGS|^ALOX|^PTGES|^HPGDS$|^TBXAS1$|^NOS[0-9]|^CA[0-9]|"
                         "^HSD|^SRD5A|^STS$' - note NAMPT is deliberately absent (NAD synthesis is "
                         'biosynthetic) and transporters (SLC*, ABC*, ATP1A/B*, AQP*) are deliberately absent '
                         'too, because SLC7A11 and SLC1A5 genuinely limit precursor supply.')
    ap.add_argument('--exclude-target-genes', default='',
                    help='drop a drug from the matched group if its matching target hits this pattern. '
                         'Needed because GEM subsystem membership is bookkeeping, not mechanism: of the 44 '
                         'PRISM drugs whose inhibited target sits in a nucleotide subsystem, 9 are '
                         'PHOSPHODIESTERASE inhibitors (papaverine, anagrelide, enoximone, cilostamide, '
                         'zardaverine, YM-976, trequinsin, avanafil, arctigenin) - cAMP/cGMP hydrolysis is '
                         'purine metabolism in the GEM, but a PDE inhibitor does not deplete nucleotide pools. '
                         'Others in the same category: GCK/HK1 (hexokinase, glycolytic), EEF2 (GTP-dependent '
                         'translation), TXNRD1/2 (thioredoxin reductase), XDH. A sensible value is '
                         "'^PDE|^TXNRD|^XDH|^GCK$|^HK1$|^EEF2$'.")
    ap.add_argument('--target-effect', default='inh',
                    help="required effect on the matching target ('inh', 'act', or '' for any). Default "
                         "'inh': activating a biosynthetic enzyme does not predict a capacity collapse.")
    ap.add_argument('--group-jaccard', type=float, default=0.0,
                    help='if > 0, collapse drugs into mechanism groups by Jaccard similarity of their target '
                         'vectors (single-linkage at this threshold) and correlate GROUP MEDIANS instead of '
                         'individual drugs. This fixes pseudo-replication, which is real here: among 747 '
                         'PRISM drugs with annotated targets there are only 593 distinct target vectors, and '
                         '237 drugs share a vector with another - 8 MEK inhibitors as {MAP2K1,MAP2K2}, 7 '
                         'EGFR, 6 TOP2A, 6 PI3K, 5 tubulin, 5 HSP90, 5 mTOR. Treating those as independent '
                         'observations inflates every p-value reported without it.')
    ap.add_argument('--group-report', type=int, default=0,
                    help='if > 0, report the differential WITHIN each Jaccard mechanism group that has at '
                         'least this many drugs, for the set named by --group-report-set. This is the '
                         'stratified counterpart to --group-jaccard: instead of collapsing groups to kill '
                         'pseudo-replication, it asks WHICH mechanisms carry an association. Needs --targets. '
                         'The large groups are mostly non-metabolic (MET 39 drugs, ATR/MTOR 35, KIT/PDGFR 24, '
                         'HDAC 16, MEK 14, PIM1 14, NR3C1 12 at Jaccard 0.5), so they double as a '
                         'negative-control panel: under the hypothesis they should show nothing.')
    ap.add_argument('--group-report-set', default='dntp_endpoints')
    ap.add_argument('--group-report-variant', default='all')
    ap.add_argument('--group-report-jaccard', type=float, default=0.5)
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

    # Target-based selection: principled where the MoA regex is free text.
    sel_drugs, sel_label, tvectors = None, '', {}
    if args.targets:
        if not args.nodes:
            raise SystemExit('--targets needs --nodes (data/networks/combined/nodes.tsv) for the '
                             'gene -> GEM subsystem map')
        g2s = gem_gene_subsystems(args.nodes)
        dte = drug_target_effects(args.targets)
        globals()['dte'] = dte
        tvectors = {d: frozenset(dte.get(d.lower(), {})) for d in pheno}
        sub_re = re.compile(args.target_subsystems, re.I) if args.target_subsystems else None
        exc_re = re.compile(args.exclude_target_genes) if args.exclude_target_genes else None
        cls_re = re.compile(args.exclude_target_classes) if args.exclude_target_classes else None
        srcmax = {}
        if args.min_target_sources > 1:
            with open(args.targets, newline='') as fh:
                for r in csv.DictReader(fh, delimiter='\t'):
                    n = (r.get('pertname') or '').strip().lower()
                    if not n:
                        continue
                    try:
                        k = int(r.get('n_sources') or 0)
                    except ValueError:
                        k = 0
                    key = (n, r['target'])
                    srcmax[key] = max(srcmax.get(key, 0), k)
        relev = {}
        if args.relevance:
            with open(args.relevance, newline='') as fh:
                for r in csv.DictReader(fh, delimiter='\t'):
                    relev.setdefault(r['set'], {})[r['gene']] = int(r['hops'])
            print('relevance table: %d sets, max hops allowed %d' % (len(relev), args.max_hops))
        sel_drugs, n_excluded = set(), 0
        for d in pheno:
            hit = False
            for gene, eff in dte.get(d.lower(), {}).items():
                subs = g2s.get(gene)
                if not subs:
                    continue
                if sub_re is not None and not any(sub_re.search(x) for x in subs):
                    continue
                if args.target_effect and eff != args.target_effect:
                    continue
                if exc_re is not None and exc_re.search(gene):
                    continue
                if cls_re is not None and cls_re.search(gene):
                    continue
                if srcmax and srcmax.get((d.lower(), gene), 0) < args.min_target_sources:
                    continue
                hit = True
                break
            if hit:
                sel_drugs.add(d)
            elif exc_re is not None:
                # would it have matched without the exclusion?
                for gene, eff in dte.get(d.lower(), {}).items():
                    subs = g2s.get(gene)
                    if subs and (sub_re is None or any(sub_re.search(x) for x in subs)) \
                            and (not args.target_effect or eff == args.target_effect):
                        n_excluded += 1
                        break
        if exc_re is not None:
            print('excluded %d drugs whose only matching target hits %r (GEM subsystem membership is '
                  'bookkeeping, not mechanism)' % (n_excluded, args.exclude_target_genes))
        sel_label = 'target in a matching GEM subsystem%s' % (
            ', effect=%s' % args.target_effect if args.target_effect else '')
        print('annotated targets: %d of %d drugs have any; %d have a GEM target in a subsystem matching '
              'the pattern%s -> these are the MATCHING group, the rest are the control'
              % (sum(1 for d in pheno if tvectors.get(d)), len(pheno), len(sel_drugs),
                 ' with effect=%s' % args.target_effect if args.target_effect else ''))
    if args.group_jaccard > 0 and not tvectors:
        if not args.targets:
            raise SystemExit('--group-jaccard needs --targets to build the target vectors')
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

    expr_scores, expr_map = {}, {}
    if args.expr_run:
        for kv in args.expr_set_map.split(','):
            if '=' in kv:
                a, b = kv.split('=', 1)
                expr_map[a.strip()] = b.strip()
        for c in list(scores):
            q = os.path.join(args.scores_root, c, args.expr_run, 'per_drug_scores.tsv')
            if os.path.isfile(q):
                expr_scores[c] = read_scores(q)
        if len(expr_scores) < 2:
            print('expression run %r found for %d cell line(s) - the conditional test needs at least 2, '
                  'so it is skipped. Run the expr arm: arm=expr bash scripts/mc_phenotype_assoc/'
                  'launch_prism_cell.sh' % (args.expr_run, len(expr_scores)))
            expr_scores = {}
        else:
            print('expression scores found for: %s  -> conditional test enabled'
                  % ', '.join(sorted(expr_scores)))

    def _expr_key(set_name, statistic):
        return (expr_map.get(set_name, set_name), args.expr_variant, statistic)

    def _partial(a, b, c):
        """Partial Spearman of a vs b controlling for c, on ranks; returns (r, p)."""
        ra, rb, rc = (sps.rankdata(x) for x in (a, b, c))
        def pear(x, y):
            x = x - x.mean(); y = y - y.mean()
            d = np.sqrt((x * x).sum() * (y * y).sum())
            return float((x * y).sum() / d) if d > 0 else np.nan
        rab, rac, rbc = pear(ra, rb), pear(ra, rc), pear(rb, rc)
        den = np.sqrt((1 - rac ** 2) * (1 - rbc ** 2))
        if not np.isfinite(den) or den <= 0:
            return np.nan, np.nan
        r = (rab - rac * rbc) / den
        n = len(a)
        if n <= 4 or not np.isfinite(r) or abs(r) >= 1:
            return r, np.nan
        t = r * np.sqrt((n - 3) / max(1e-12, 1 - r ** 2))
        return r, float(2 * sps.t.sf(abs(t), n - 3))

    def _is_match(d):
        if sel_drugs is not None:
            return d in sel_drugs
        if moa_re is not None:
            return bool(moa_re.search(moa.get(d, '')))
        return False

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
            # Mechanism grouping: collapse each Jaccard group to its median before
            # correlating, so 8 MEK inhibitors count once rather than eight times.
            # A group counts as MATCHING if any of its members does.
            if args.group_jaccard > 0:
                idx = {d: i for i, d in enumerate(drugs)}
                gdP, gdS, gmatch, gnames = [], [], [], []
                for grp in jaccard_groups({d: tvectors.get(d, frozenset()) for d in drugs},
                                          args.group_jaccard):
                    ii = [idx[d] for d in grp]
                    gdP.append(float(np.median(dP[ii])))
                    gdS.append(float(np.median(dS[ii])))
                    gmatch.append(any(_is_match(d) for d in grp))
                    gnames.append(grp[0] if len(grp) == 1 else '%s+%d' % (grp[0], len(grp) - 1))
                dP, dS = np.array(gdP), np.array(gdS)
                drugs, group_match = gnames, np.array(gmatch)
            else:
                group_match = None
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
            if group_match is None:
                naive = corr(np.array([sx[d] for d in drugs]),
                             np.array([pheno[d][X] for d in drugs]))[0]
            else:
                naive = np.nan   # undefined once drugs are collapsed into groups

            r_moa = p_moa = r_oth = p_oth = np.nan
            n_moa = n_oth = 0
            m = None
            if group_match is not None:
                m = group_match
            elif sel_drugs is not None:
                # --relevance narrows the matched group PER SET: a drug counts only if one
                # of its (already filtered) inhibited targets lies within --max-hops of a
                # producing route to THIS set. Hence the mask is rebuilt for each set.
                if relev:
                    hops = relev.get(key[0], {})
                    ok_rel = set()
                    for d in sel_drugs:
                        for gene in dte.get(d.lower(), {}):
                            h = hops.get(gene)
                            if h is not None and h <= args.max_hops:
                                ok_rel.add(d)
                                break
                    m = np.array([d in ok_rel for d in drugs])
                else:
                    m = np.array([d in sel_drugs for d in drugs])
            elif moa_re is not None:
                m = np.array([bool(moa_re.search(moa.get(d, ''))) for d in drugs])

            # Conditional on expression: the differential of the matching GENE set, and
            # the MetChange differential with that gene set partialled out.
            r_expr = p_expr = r_part = p_part = np.nan
            r_expr_moa = r_part_moa = p_part_moa = np.nan
            r_expr_oth = r_part_oth = np.nan
            n_cond = n_cond_moa = 0
            if expr_scores and group_match is None and X in expr_scores and Y in expr_scores:
                ek = _expr_key(key[0], key[2])
                ex, ey = expr_scores[X].get(ek), expr_scores[Y].get(ek)
                if ex and ey:
                    cd = [d for d in drugs if d in ex and d in ey
                          and np.isfinite(ex[d]) and np.isfinite(ey[d])]
                    if len(cd) >= args.min_drugs:
                        n_cond = len(cd)
                        ii = {d: i for i, d in enumerate(drugs)}
                        sub = [ii[d] for d in cd]
                        eP, eS = dP[sub], dS[sub]
                        eE = np.array([ex[d] - ey[d] for d in cd])
                        if np.nanstd(eE) > 0:
                            r_expr, p_expr = corr(eE, eP)
                            r_part, p_part = _partial(eS, eP, eE)
                    # The pre-registered statistic is the partial WITHIN the matched
                    # subgroup, not over all drugs: the all-drug partial answers a
                    # different question (does expression explain the pooled signal,
                    # which is itself confounded and near-null).
                    if m is not None:
                        mm = [d for d in cd if _is_match(d)]
                        oo = [d for d in cd if not _is_match(d)]
                        for lst, tag in ((mm, 'moa'), (oo, 'oth')):
                            if len(lst) < max(8, args.moa_min_drugs // 2):
                                continue
                            jj = [ii[d] for d in lst]
                            aS, aP = dS[jj], dP[jj]
                            aE = np.array([ex[d] - ey[d] for d in lst])
                            if np.nanstd(aE) <= 0 or np.nanstd(aS) <= 0:
                                continue
                            re_, _ = corr(aE, aP)
                            rp_, pp_ = _partial(aS, aP, aE)
                            if tag == 'moa':
                                r_expr_moa, r_part_moa, p_part_moa, n_cond_moa = re_, rp_, pp_, len(lst)
                            else:
                                r_expr_oth, r_part_oth = re_, rp_
            # MoA strata: the matching group is where the hypothesis applies at all,
            # the complement is the built-in control.
            if m is not None:
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
                         'n_other': n_oth, 'r_other': r_oth, 'p_other': p_oth,
                         'n_conditional': n_cond, 'r_expr': r_expr, 'p_expr': p_expr,
                         'r_partial_mc_given_expr': r_part, 'p_partial': p_part,
                         'n_cond_moa': n_cond_moa, 'r_expr_moa': r_expr_moa,
                         'r_partial_moa': r_part_moa, 'p_partial_moa': p_part_moa,
                         'r_expr_other': r_expr_oth, 'r_partial_other': r_part_oth})
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
            print('STRATIFIED by %s: matching (where the hypothesis applies) vs complement (control)'
                  % (sel_label or ('MoA %r' % args.moa_regex)))
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
    if expr_scores:
        ok = [r for r in rows if np.isfinite(r['r_partial_mc_given_expr'])]
        if ok:
            print()
            print('CONDITIONAL ON EXPRESSION: does the MetChange differential survive controlling for the '
                  'matching gene set?')
            h = ('%-9s %-20s %-11s %5s %9s %9s %11s %9s'
                 % ('pair', 'set', 'variant', 'n', 'r_MetChg', 'r_expr', 'r_MC|expr', 'p'))
            print(h)
            print('-' * len(h))
            for r in sorted(ok, key=lambda r: r['p_partial'])[:args.top]:
                print('%-9s %-20s %-11s %5d %+9.3f %+9.3f %+11.3f %9.3g'
                      % ('%s/%s' % (r['cell_x'], r['cell_y']), r['set'][:20], r['variant'],
                         r['n_conditional'], r['r_diff'], r['r_expr'],
                         r['r_partial_mc_given_expr'], r['p_partial']))
            rm = np.array([r['r_diff'] for r in ok])
            re_ = np.array([r['r_expr'] for r in ok])
            rp = np.array([r['r_partial_mc_given_expr'] for r in ok])
            print()
            print('over %d combinations: median |r_MetChange| %.3f, |r_expression| %.3f, '
                  '|r_MetChange|expression| %.3f'
                  % (len(ok), np.nanmedian(np.abs(rm)), np.nanmedian(np.abs(re_)),
                     np.nanmedian(np.abs(rp))))
            keep = np.nanmedian(np.abs(rp)) / max(1e-9, np.nanmedian(np.abs(rm)))
            print('the MetChange association retains %.0f%% of its magnitude once expression is partialled '
                  'out (100%% = expression explains none of it, 0%% = all of it)' % (100 * keep))
    if args.group_report:
        if not tvectors:
            print('\n--group-report needs --targets')
        else:
            key = (args.group_report_set, args.group_report_variant, args.statistic)
            avail = [c for c in scores if key in scores[c]]
            if len(avail) < 2:
                print('\n--group-report: %s not scored in >=2 cells' % (key,))
            else:
                print()
                print('WITHIN-GROUP DIFFERENTIAL for %s %s, Jaccard >= %.2f, groups with >= %d drugs'
                      % (args.group_report_set, args.group_report_variant,
                         args.group_report_jaccard, args.group_report))
                print('Each row: the correlation computed ACROSS that group\'s drugs, per cell pair, then')
                print('summarised over pairs. "matched" = the group targets a selected GEM subsystem.')
                base = sorted(set(pheno) & set.intersection(*[set(scores[c][key]) for c in avail]))
                grps = [g for g in jaccard_groups({d: tvectors.get(d, frozenset()) for d in base},
                                                  args.group_report_jaccard)
                        if len(g) >= args.group_report]
                h = ('%-30s %5s %4s %9s %9s %8s'
                     % ('group (shared targets)', 'n', 'mtch', 'median r', 'frac r>0', 'npairs'))
                print(h)
                print('-' * len(h))
                out = []
                for g in sorted(grps, key=len, reverse=True):
                    rs = []
                    for X, Y in itertools.combinations(sorted(avail), 2):
                        if X not in pheno[g[0]] if g[0] in pheno else True:
                            pass
                        dd = [d for d in g if X in pheno.get(d, {}) and Y in pheno.get(d, {})
                              and d in scores[X][key] and d in scores[Y][key]
                              and np.isfinite(scores[X][key][d]) and np.isfinite(scores[Y][key][d])]
                        if len(dd) < max(8, args.group_report // 2):
                            continue
                        a = np.array([scores[X][key][d] - scores[Y][key][d] for d in dd])
                        b = np.array([pheno[d][X] - pheno[d][Y] for d in dd])
                        if np.nanstd(a) == 0 or np.nanstd(b) == 0:
                            continue
                        rs.append(corr(a, b)[0])
                    if len(rs) < 3:
                        continue
                    rs = np.array([x for x in rs if np.isfinite(x)])
                    tv = sorted(tvectors.get(g[0], frozenset()))
                    lbl = ','.join(tv[:3]) + ('...' if len(tv) > 3 else '') or '(no targets)'
                    mtch = 'yes' if any(_is_match(d) for d in g) else 'no'
                    out.append((len(g), lbl, mtch, float(np.median(rs)), float(np.mean(rs > 0)), len(rs)))
                for n_, lbl, mtch, med, frac, npair in out:
                    print('%-30s %5d %4s %+9.3f %8.0f%% %8d' % (lbl[:30], n_, mtch, med, 100 * frac, npair))
                if out:
                    mm = [o for o in out if o[2] == 'yes']
                    nn = [o for o in out if o[2] == 'no']
                    print()
                    print('median of the group medians: matched groups %+.3f (n=%d), '
                          'non-matched %+.3f (n=%d)'
                          % (np.median([o[3] for o in mm]) if mm else np.nan, len(mm),
                             np.median([o[3] for o in nn]) if nn else np.nan, len(nn)))
                    print('NOTE the non-matched groups are mostly kinase/HDAC/receptor mechanisms and act as '
                          'a negative-control panel here.')
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
