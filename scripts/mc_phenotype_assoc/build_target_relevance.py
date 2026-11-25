#!/usr/bin/env python3
"""Which drug targets can actually perturb a metabolite set's production?

Emits `gene, set, hops`: the shortest number of reaction steps from a gene's
reactions, following PRODUCTION only, to a member of the named metabolite set.

WHY THIS EXISTS. Selecting "metabolically-acting" drugs by the GEM **subsystem** of
their target does not work, and has failed three separate times in report 10 by
producing a plausible-looking drug group that was mostly irrelevant:

  - PI3K qualified as metabolic because phosphoinositides are metabolites, so 11
    PI3K inhibitors entered the group (report 10, 10.1);
  - nine PHOSPHODIESTERASE inhibitors qualified under purine metabolism, because
    cAMP/cGMP hydrolysis is purine metabolism (11.3);
  - thirteen drugs qualified under amino-acid metabolism via **CYP2C19**, which
    metabolises them rather than being their site of action, plus six MAO
    inhibitors (monoamine catabolism sits in tyrosine/tryptophan) and two histone
    methyltransferase inhibitors (their targets consume SAM, so they live in
    cysteine/methionine metabolism) (13.4).

Human-GEM is not wrong in any of these cases. The subsystem a gene's reactions
belong to simply does not say whether a drug acting on that gene perturbs flux
toward a biomass precursor. **Direction does.** Requiring the target to lie on a
route that PRODUCES a set member excludes PDE, PI3K and EZH2 automatically - they
are never reached - which is exactly the hand-made exclusion list of 11.3,
derived from the stoichiometry instead of by inspection.

WHAT IT CANNOT DO, measured rather than assumed. Hop distance does **not** separate
biosynthetic routes from incidental ones. Distances to a dNTP:

    NME1/2, RRM1/2            1        the NDP kinases and ribonucleotide reductase
    TK1, TYMS, EHMT2          4        ...and a histone methyltransferase
    DHFR, IMPDH1, TYMP, ADA,
      CYP2C19                 5        ...and a drug-metabolising CYP
    DHODH, MAOA               6        ...and monoamine oxidase
    PDE*, PIK3*, EZH2         never reached

So any cutoff deep enough to include DHODH also includes MAO-A; any that includes
DHFR also includes CYP2C19; any that includes TYMS also includes G9a. The gene
count explodes with depth too - 19 genes at hop 1, 402 at 4, 837 at 5, 1,415 at 6
of 2,848. **Use this filter together with the class exclusions and the annotation
confidence threshold documented in report 10, 14, not instead of them.**

Currency metabolites are excluded from the traversal, for the same reason the
network builder removes them: ATP and H2O appear in thousands of reactions, so a
path through them reaches everything in two steps and the distance stops meaning
anything.
"""
import argparse, collections, csv, os, re, sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Mirrors CURRENCY_NAMES in build_metabolite_sets.py, plus the nucleoside
# tri/diphosphates, which otherwise short-circuit every nucleotide route.
CURRENCY = {'ATP', 'ADP', 'AMP', 'GTP', 'GDP', 'GMP', 'UTP', 'UDP', 'UMP', 'CTP', 'CDP', 'CMP',
            'Pi', 'PPi', 'H+', 'H2O', 'NAD+', 'NADH', 'NADP+', 'NADPH', 'CO2', 'O2', 'CoA',
            'NH3', 'NH4+', 'Na+', 'K+', 'Cl-', 'H2O2', 'FAD', 'FADH2', 'FMN', 'SAM', 'SAH'}


def parse_gem(path, genes_tsv):
    """-> (metabolite id -> name, reaction id -> (coeffs, [gene symbols]))"""
    ens2sym = {}
    if os.path.exists(genes_tsv):
        with open(genes_tsv) as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                g = r.get('genes')
                s = r.get('geneSymbols') or ''
                if g and s:
                    ens2sym[g] = s
    mets, rxns = {}, {}
    re_id = re.compile(r'^\s+- id: "(.*)"')
    re_nm = re.compile(r'^\s+- name: "(.*)"')
    re_met = re.compile(r'^\s+- (MAM\w+): (-?[\d.eE+]+)')
    re_gpr = re.compile(r'^\s+- gene_reaction_rule: "(.*)"')
    sect = cur = rid = None
    coeffs, gpr = {}, ''
    with open(path) as fh:
        for line in fh:
            if line.startswith('- metabolites:'):
                sect = 'm'; continue
            if line.startswith('- reactions:'):
                sect = 'r'; continue
            if line.startswith('- genes:'):
                break
            if sect == 'm':
                m = re_id.match(line)
                if m:
                    cur = m.group(1); continue
                m = re_nm.match(line)
                if m and cur:
                    mets[cur] = m.group(1)
            elif sect == 'r':
                m = re_id.match(line)
                if m:
                    if rid:
                        rxns[rid] = (coeffs, gpr)
                    rid, coeffs, gpr = m.group(1), {}, ''
                    continue
                m = re_met.match(line)
                if m:
                    coeffs[m.group(1)] = float(m.group(2)); continue
                m = re_gpr.match(line)
                if m:
                    gpr = m.group(1)
    if rid:
        rxns[rid] = (coeffs, gpr)
    out = {}
    for r, (c, g) in rxns.items():
        out[r] = (c, [ens2sym.get(x, x) for x in re.findall(r'ENSG\d+', g)])
    print('GEM: %d metabolites, %d reactions, %d ENSG->symbol' % (len(mets), len(out), len(ens2sym)),
          file=sys.stderr)
    return mets, out


def read_sets(path):
    """set name -> {metabolite names}, from the `all` variant (compartment-agnostic)."""
    out = collections.defaultdict(set)
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            if r['variant'] != 'all':
                continue
            for mid in r['members'].split(','):
                if mid.strip():
                    out[r['set']].add(mid.strip())
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--gem', default=os.path.join(REPO, 'data/GEMs/Human-GEM/model/Human-GEM.yml'))
    ap.add_argument('--genes', default=os.path.join(REPO, 'data/GEMs/Human-GEM/model/genes.tsv'))
    ap.add_argument('--sets', default=os.path.join(REPO, 'results/mc_phenotype_assoc/metabolite_sets/metabolite_sets.tsv'))
    ap.add_argument('--annotation', default=os.path.join(REPO, 'results/mc_phenotype_assoc/metabolite_sets/metabolite_annotation.tsv'),
                    help='maps metabolite node ids to names; needed because the set files hold ids')
    ap.add_argument('--max-hops', type=int, default=6)
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    mets, rxns = parse_gem(args.gem, args.genes)
    name_of = dict(mets)
    sets = read_sets(args.sets)

    # the set files store node ids (MAM...c); map to bare names via the annotation
    # table when present, else assume the member strings are already names
    id2name = {}
    if os.path.exists(args.annotation):
        with open(args.annotation) as fh:
            for r in csv.DictReader(fh, delimiter='\t'):
                if r.get('metabolite') and r.get('name'):
                    id2name[r['metabolite']] = r['name']
    print('metabolite id -> name: %d' % len(id2name), file=sys.stderr)

    def member_names(members):
        return {id2name[m] for m in members if m in id2name}

    prod_by = collections.defaultdict(list)
    subs_of = collections.defaultdict(list)
    genes_of = {}
    for r, (c, g) in rxns.items():
        genes_of[r] = g
        for mid, v in c.items():
            if v > 0:
                prod_by[mid].append(r)
            elif v < 0 and name_of.get(mid, '') not in CURRENCY:
                subs_of[r].append(mid)

    by_name = collections.defaultdict(set)
    for mid, nm in name_of.items():
        by_name[nm].add(mid)

    rows, summary = [], []
    for sname, members in sorted(sets.items()):
        want = member_names(members)
        targets = {mid for nm in want for mid in by_name.get(nm, ())}
        if not targets:
            summary.append((sname, len(members), 0, 0))
            continue
        seen_g, seen_m, frontier = {}, set(targets), set(targets)
        for hop in range(1, args.max_hops + 1):
            nxt = set()
            for mid in frontier:
                for r in prod_by.get(mid, []):
                    for gene in genes_of.get(r, []):
                        seen_g.setdefault(gene, hop)
                    for s in subs_of.get(r, []):
                        if s not in seen_m:
                            seen_m.add(s); nxt.add(s)
            frontier = nxt
            if not frontier:
                break
        for gene, hop in sorted(seen_g.items()):
            rows.append({'gene': gene, 'set': sname, 'hops': hop})
        summary.append((sname, len(members), len(targets), len(seen_g)))

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or '.', exist_ok=True)
    with open(args.out, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=['gene', 'set', 'hops'], delimiter='\t')
        w.writeheader()
        w.writerows(rows)
    print('%-24s %8s %8s %8s' % ('set', 'members', 'matched', 'genes'), file=sys.stderr)
    for sname, nmem, ntgt, ngene in summary:
        print('%-24s %8d %8d %8d' % (sname, nmem, ntgt, ngene), file=sys.stderr)
    print('\nwrote %s (%d rows)' % (args.out, len(rows)), file=sys.stderr)


if __name__ == '__main__':
    main()
