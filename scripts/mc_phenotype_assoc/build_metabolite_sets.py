#!/usr/bin/env python3
"""Build metabolite sets for the MetChange / mitochondrial-phenotype association.

Sets come from Human-GEM subsystems: a set is every metabolite taking part in a
reaction of that subsystem. Compartment variants are emitted because the
MetChange features are compartment-specific (MAM00003c vs MAM00003m), and
respiration is a mitochondrial phenotype while glycolysis is cytosolic -
restricting the compartment is what makes the sets specific rather than just
"central metabolism".

Hub metabolites (ATP, H+, H2O, Pi, NAD+ ...) take part in hundreds of reactions,
so a subsystem set that contains them is partly a set of the hubs, and a
difference in it can reflect global network capacity rather than that pathway.
Every set is therefore emitted twice, with and without hubs (`_nohub`), and the
hub list is written out so the choice is visible.

Also emits curated endpoint sets. For a respiration phenotype the informative
features are the endpoints of the pathway (O2 and NADH consumed, lactate
produced), not every intermediate: at higher flux the intermediates of a pathway
often stay flat while throughput rises (report 06, section 7).

Note the currency metabolites (ATP, NAD+/NADH, O2, ...) ARE MetChange features
but are NOT network nodes, so they can be used in this label-level test and
cannot be used to build a model query (report 07, section 5).

`--emit genes` instead writes the GPR **gene** sets of the same subsystems, taken
from each reaction's `gene_reaction_rule`, so the identical association test can
be run on the expression labels as the baseline that says whether the
stoichiometric reading adds anything over reading the genes. The GEM states those
rules in Ensembl IDs, which is also what the expression label matrices are indexed
by, so no identifier mapping is involved. Genes have no compartment, so the
variants there are `gene` and `gene_nohub` (hub = a gene used by very many
reactions, the analogue of a currency metabolite).

Output: `metabolite_sets.tsv` (or `gene_sets.tsv`) with columns
set / variant / n_members / members, plus `metabolite_annotation.tsv` /
`gene_annotation.tsv` (id, name, compartment, GEM degree, is_hub) and
`sets_summary.txt`.
"""
import argparse, collections, os, re, sys

# Human-GEM subsystem names, exactly as they appear in the model
SUBSYSTEMS = {
    'oxphos':     ['Oxidative phosphorylation'],
    'tca':        ['Tricarboxylic acid cycle and glyoxylate/dicarboxylate metabolism'],
    'glycolysis': ['Glycolysis / Gluconeogenesis'],
    'pyruvate':   ['Pyruvate metabolism'],
    'ppp':        ['Pentose phosphate pathway', 'Pentose Phosphate Pathway'],
}

# Curated endpoint sets, by metabolite name as Human-GEM spells it.
CURATED = {
    # what an ETC block / uncoupling acts on directly
    'etc_endpoints':        ['O2', 'NADH', 'NAD+', 'FADH2', 'FAD', 'ubiquinol', 'ubiquinone',
                             'H2O', 'H+'],
    # the energy charge: for the membrane-potential phenotype (Attene-Ramos)
    'energy_currency':      ['ATP', 'ADP', 'AMP', 'NAD+', 'NADH', 'NADP+', 'NADPH',
                             'Pi', 'O2', 'CoA'],
    'adenylates':           ['ATP', 'ADP', 'AMP', 'Pi'],
    'redox_couples':        ['NAD+', 'NADH', 'NADP+', 'NADPH'],
    # throughput endpoints rather than intermediates
    'glycolysis_endpoints': ['glucose', 'pyruvate', 'L-lactate'],
    'tca_endpoints':        ['acetyl-CoA', 'citrate', 'AKG', 'succinate',
                             'fumarate', 'malate', 'CO2'],
}

# Hubs by name, in addition to whatever the degree cutoff picks up
CURRENCY_NAMES = {'ATP', 'ADP', 'AMP', 'NAD+', 'NADH', 'NADP+', 'NADPH', 'H+', 'H2O',
                  'Pi', 'PPi', 'CoA', 'O2', 'CO2', 'NH3', 'HCO3-', 'FAD', 'FADH2',
                  'GTP', 'GDP', 'UTP', 'UDP', 'S-adenosyl-L-methionine'}


def parse_gem(path):
    """id -> (name, compartment), and reaction -> (subsystem, [metabolite ids]).

    A line parser, not a YAML load: the file is 7.7 MB of regular !!omap blocks
    and PyYAML's omap handling is both slow and fussy here.
    """
    mets, rxns = {}, []
    re_gpr = re.compile(r'^\s+- gene_reaction_rule: "(.*)"')
    section = None
    cur, cur_mets, in_mets = {}, [], False
    re_id = re.compile(r'^\s+- id: "(.*)"')
    re_name = re.compile(r'^\s+- name: "(.*)"')
    re_comp = re.compile(r'^\s+- compartment: "(.*)"')
    re_sub = re.compile(r'^\s+- subsystem: "(.*)"')
    re_met = re.compile(r'^\s+- (MAM\w+): (-?[\d.eE+]+)')

    def flush():
        if not cur:
            return
        if section == 'metabolites' and 'id' in cur:
            mets[cur['id']] = (cur.get('name', ''), cur.get('compartment', ''))
        elif section == 'reactions' and 'id' in cur:
            rxns.append((cur['id'], cur.get('subsystem', ''), list(cur_mets),
                         re.findall(r'ENSG\d+', cur.get('gpr', ''))))

    with open(path) as fh:
        for line in fh:
            if line.startswith('- metabolites:'):
                flush(); section, cur, cur_mets = 'metabolites', {}, []
                continue
            if line.startswith('- reactions:'):
                flush(); section, cur, cur_mets = 'reactions', {}, []
                continue
            if line.startswith('- genes:') or line.startswith('- compartments:'):
                flush(); section, cur, cur_mets = None, {}, []
                continue
            if section is None:
                continue
            if line.strip() == '- !!omap':
                flush(); cur, cur_mets, in_mets = {}, [], False
                continue
            if re.match(r'^\s+- metabolites: !!omap', line):
                in_mets = True
                continue
            m = re_met.match(line)
            if m and in_mets:
                cur_mets.append(m.group(1))
                continue
            for key, rx in (('id', re_id), ('name', re_name),
                            ('compartment', re_comp), ('subsystem', re_sub),
                            ('gpr', re_gpr)):
                m2 = rx.match(line)
                if m2:
                    cur[key] = m2.group(1)
                    if key != 'id':
                        in_mets = False
                    break
    flush()
    return mets, rxns


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap.add_argument('--gem', default=os.path.join(repo, 'data/GEMs/Human-GEM/model/Human-GEM.yml'))
    ap.add_argument('--features', help='optional: a label matrix (log2FC.csv) whose first '
                                      'column lists the MetChange features; sets are '
                                      'restricted to features that actually exist')
    ap.add_argument('--hub-percentile', type=float, default=99.0,
                    help='metabolites above this percentile of GEM reaction degree are hubs')
    ap.add_argument('--min-set-size', type=int, default=5)
    ap.add_argument('--emit', choices=['metabolites', 'genes'], default='metabolites',
                    help='genes = the subsystems\' GPR genes (Ensembl), for the expression baseline')
    ap.add_argument('--out-dir', default=os.path.join(repo, 'data/validation/metabolite_sets'))
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    mets, rxns = parse_gem(args.gem)
    print('GEM: %d metabolites, %d reactions' % (len(mets), len(rxns)), file=sys.stderr)

    present = None
    if args.features:
        present = set()
        with open(args.features) as fh:
            next(fh)
            for line in fh:
                present.add(line.split(',', 1)[0].strip().strip('"'))
        print('label matrix: %d features' % len(present), file=sys.stderr)

    # degree over the whole model, for the hub cutoff
    degree = collections.Counter()
    for _, _, ms, gs in rxns:
        for m in set(gs if args.emit == 'genes' else ms):
            degree[m] += 1
    import numpy as np
    cutoff = float(np.percentile(list(degree.values()), args.hub_percentile)) if degree else 0.0
    if args.emit == 'genes':
        universe = sorted(degree)
        mets = {g: (g, '') for g in universe}
        is_hub = {g: degree[g] >= cutoff for g in universe}
    else:
        is_hub = {m: (degree[m] >= cutoff or mets.get(m, ('', ''))[0] in CURRENCY_NAMES)
                  for m in mets}
    print('hub cutoff: degree >= %.0f (%d hubs of %d)'
          % (cutoff, sum(is_hub.values()), len(is_hub)), file=sys.stderr)

    by_sub = collections.defaultdict(set)
    for _, sub, ms, gs in rxns:
        if sub:
            by_sub[sub].update(gs if args.emit == 'genes' else ms)

    raw = {}
    for key, names in SUBSYSTEMS.items():
        members = set()
        for n in names:
            members |= by_sub.get(n, set())
        if not members:
            print('WARNING: subsystem %r matched nothing' % names, file=sys.stderr)
        raw[key] = members
    name_to_ids = collections.defaultdict(set)
    for mid, (nm, _) in mets.items():
        name_to_ids[nm].add(mid)
    for key, names in ({} if args.emit == 'genes' else CURATED).items():
        members = set()
        for n in names:
            if n not in name_to_ids:
                print('WARNING: curated metabolite %r not in the GEM' % n, file=sys.stderr)
            members |= name_to_ids.get(n, set())
        raw[key] = members

    rows, summary = [], []
    for key, members in sorted(raw.items()):
        comps = ((('gene', None),) if args.emit == 'genes'
                 else (('all', None), ('mito', 'm'), ('cyto', 'c')))
        for comp_label, comp in comps:
            sel = {m for m in members
                   if comp is None or mets.get(m, ('', ''))[1] == comp}
            for hub_label, keep in (('', True), ('_nohub', False)):
                s = sel if keep else {m for m in sel if not is_hub.get(m)}
                if present is not None:
                    s = {m for m in s if m in present}
                if len(s) < args.min_set_size:
                    continue
                variant = comp_label + hub_label
                rows.append((key, variant, len(s), ','.join(sorted(s))))
                summary.append('%-22s %-10s %4d' % (key, variant, len(s)))

    out = os.path.join(args.out_dir,
                       'gene_sets.tsv' if args.emit == 'genes' else 'metabolite_sets.tsv')
    with open(out, 'w') as fh:
        fh.write('set\tvariant\tn_members\tmembers\n')
        for r in rows:
            fh.write('%s\t%s\t%d\t%s\n' % r)
    ann_path = os.path.join(args.out_dir,
                            'gene_annotation.tsv' if args.emit == 'genes'
                            else 'metabolite_annotation.tsv')
    with open(ann_path, 'w') as fh:
        fh.write('%s\tname\tcompartment\tgem_degree\tis_hub\tin_labels\n'
                 % ('gene' if args.emit == 'genes' else 'metabolite'))
        for mid, (nm, comp) in sorted(mets.items()):
            fh.write('%s\t%s\t%s\t%d\t%s\t%s\n'
                     % (mid, nm, comp, degree[mid], int(bool(is_hub.get(mid))),
                        '' if present is None else int(mid in present)))
    with open(os.path.join(args.out_dir, 'sets_summary%s.txt'
                           % ('_genes' if args.emit == 'genes' else '')), 'w') as fh:
        fh.write('set                    variant    size\n')
        fh.write('\n'.join(summary) + '\n')
    print('\n'.join(summary))
    print('wrote %s (%d sets)' % (out, len(rows)), file=sys.stderr)


if __name__ == '__main__':
    main()
