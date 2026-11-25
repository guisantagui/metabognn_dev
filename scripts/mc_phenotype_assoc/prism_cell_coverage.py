#!/usr/bin/env python3
"""Which cell lines are worth a MetChange label run? Reproduces report 10, section 15.

Answers, before any label run is launched: for each cell line present in both the
LINCS per-cell perturbation table and PRISM, how many of the drugs a given finding
rests on were actually profiled there?

WHY IT EXISTS. Report 10, 12.3 claimed "93 PRISM cell lines are in LINCS without
labels", taken from `cellinfo_beta.txt` membership. The real figure is **two**. The
non-core LINCS panels are modally 120 or 723 perturbations - a small fixed plate
from a different chemical space - and carry 5 or fewer of the 31 nucleotide-targeting
drugs. U2OS is the instructive case: 20,836 perturbations and still only 11 of the
31, so panel SIZE is not the criterion either; composition is.

This was the fifth instance of the same error in the project (naia 683 -> 10, xu
53 -> 11, the PRISM subsystem comparators, Hallinger). The rule it enforces: count
the units that can actually carry the prediction, for the specific comparison
proposed, before quoting an opportunity.
"""
import argparse, collections, csv, os, re, sys

csv.field_size_limit(10 ** 9)
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--bycell', default=os.path.join(REPO, 'data/LINCS/lincs_unique_perts_bycell.tsv'))
    ap.add_argument('--perts', default=os.path.join(REPO, 'data/LINCS/lincs_unique_perts.tsv'))
    ap.add_argument('--prism', default=os.path.join(
        REPO, 'data/validation/prism/secondary-screen-dose-response-curve-parameters.csv'))
    ap.add_argument('--phenotype', default=os.path.join(
        REPO, 'results/mc_phenotype_assoc/phenotypes_ALL/prism_viability.tsv'))
    ap.add_argument('--targets', default=os.path.join(REPO, 'data/drug_targets/combined/drug_targets.tsv'))
    ap.add_argument('--nodes', default=os.path.join(REPO, 'data/networks/combined/nodes.tsv'))
    ap.add_argument('--labels-root', default=os.path.join(REPO, 'results/labels'))
    ap.add_argument('--target-subsystems', default='nucleotide|purine|pyrimidine|folate')
    ap.add_argument('--exclude-target-classes',
                    default=r'^PDE|^PIK3|^PIP[45]K|^SPHK|^DGK|^XDH|^CYP|^UGT|^SULT|^NAT[12]$|^FMO|^EHMT'
                            r'|^EZH|^DNMT|^SETD|^KMT|^PRMT|^COMT$|^NNMT$|^MAO|^ALDH|^ADH[0-9]|^AOX1$'
                            r'|^DBH$|^SIRT|^PARP|^CD38$|^TXNRD|^GSR$|^GPX|^GST|^NQO1$|^SOD[0-9]|^CAT$'
                            r'|^PTGS|^ALOX|^PTGES|^HPGDS$|^TBXAS1$|^NOS[0-9]|^CA[0-9]|^HSD|^SRD5A|^STS$',
                    help='the class list of report 10, 14.3/14.5')
    ap.add_argument('--top', type=int, default=16)
    ap.add_argument('--out')
    args = ap.parse_args()

    from prism_differential import gem_gene_subsystems, drug_target_effects

    # LINCS: cell -> set of drug names (labels are keyed on pertname)
    lincs = collections.defaultdict(set)
    with open(args.bycell, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            lincs[r['cell'].upper()].add(r['pertname'])

    id2name = {}
    with open(args.perts, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            id2name[r['pert_id']] = r['pertname']

    brd = lambda s: '-'.join(str(s).split('-')[:2])
    prism = collections.defaultdict(set)
    with open(args.prism, newline='') as fh:
        for r in csv.DictReader(fh):
            if not r.get('auc'):
                continue
            nm = id2name.get(brd(r['broad_id']))
            if nm:
                prism[r['ccle_name'].split('_')[0].upper()].add(nm)

    # the drug set the finding rests on
    g2s = gem_gene_subsystems(args.nodes)
    dte = drug_target_effects(args.targets)
    sub = re.compile(args.target_subsystems, re.I)
    cls = re.compile(args.exclude_target_classes)
    focus = set()
    with open(args.phenotype) as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            d = r['drug']
            for gene, eff in dte.get(d.lower(), {}).items():
                if eff != 'inh' or cls.search(gene):
                    continue
                if any(sub.search(x) for x in g2s.get(gene, ())):
                    focus.add(d)
                    break
    have = {c.upper() for c in os.listdir(args.labels_root)
            if os.path.isdir(os.path.join(args.labels_root, c))}
    print('focus drug set (target in %r, effect=inh, class-excluded): %d drugs'
          % (args.target_subsystems, len(focus)))

    rows = []
    for c in sorted(set(lincs) & set(prism)):
        ov = lincs[c] & prism[c]
        rows.append({'cell': c, 'lincs_perts': len(lincs[c]), 'prism_drugs': len(prism[c]),
                     'overlap': len(ov), 'focus_covered': len(ov & focus),
                     'has_labels': c in have})
    rows.sort(key=lambda r: (-r['focus_covered'], -r['overlap']))
    print('cell lines in both tables: %d\n' % len(rows))
    print('%-14s %9s %9s %9s %12s %8s' % ('cell', 'LINCS n', 'PRISM n', 'overlap',
                                          'focus/%d' % len(focus), 'labels'))
    for r in rows[:args.top]:
        print('%-14s %9d %9d %9d %12d %8s'
              % (r['cell'], r['lincs_perts'], r['prism_drugs'], r['overlap'],
                 r['focus_covered'], 'YES' if r['has_labels'] else '-'))
    print()
    dist = collections.Counter(r['focus_covered'] for r in rows)
    print('focus-coverage distribution: %s' % dict(sorted(dist.items(), reverse=True)))
    cand = [r for r in rows if not r['has_labels'] and r['focus_covered'] >= 0.5 * len(focus)]
    print()
    if cand:
        print('WORTH A LABEL RUN (no labels yet, >= half the focus drugs):')
        for r in cand:
            print('   %-12s %d of %d focus drugs, %d PRISM drugs overlapping'
                  % (r['cell'], r['focus_covered'], len(focus), r['overlap']))
        n0 = len(have & {r['cell'] for r in rows})
        for k in range(1, len(cand) + 1):
            n = n0 + k
            print('   adding %d -> %d cells, %d pairs, Wilcoxon floor %.2g'
                  % (k, n, n * (n - 1) // 2, 2.0 ** -(n * (n - 1) // 2)))
    else:
        print('no cell line without labels covers half the focus drugs')
    if args.out:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)) or '.', exist_ok=True)
        with open(args.out, 'w', newline='') as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter='\t')
            w.writeheader(); w.writerows(rows)
        print('\nwrote %s' % args.out, file=sys.stderr)


if __name__ == '__main__':
    main()
