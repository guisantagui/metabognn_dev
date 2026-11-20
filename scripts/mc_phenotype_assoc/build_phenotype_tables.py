#!/usr/bin/env python3
"""Turn the mapped Attene-Ramos and Hallinger data into drug -> phenotype tables.

Keys on the LINCS perturbation NAME, because that is what the MetChange label
columns are named after ("<drug>_<dose>_<time>"); `fit_info_nsig.csv` carries the
parsed `drug` column that these join to.

Attene-Ramos (PubChem AID 720637, MMP in HepG2) is multiplexed: the MMP ratio
readout and a cell-viability counter screen come from the same wells, and
`Activity Summary` is the curators' combined call. A compound that kills cells
loses membrane potential trivially, so the informative contrast is MMP-down
*without* viability loss against inactive. Groups emitted:
  MMP_down   "active antagonist"                      - clean MMP loss
  cytotox    "inconclusive antagonist (cytotoxic)"    - confounded, kept separate
  MMP_up     "active agonist"                         - hyperpolarisation, its own phenotype
  inactive   "inactive"
  incon      everything else, excluded by default
`--strict-sep` additionally demands that viability is either never active or at
least `--sep-fold` times less potent than the MMP effect.

Hallinger (Dryad zkh189367 + the paper's Supplemental Table 2) is respirometry
only - all six endpoints are OCR, there is no ECAR in the deposit - so it gives
activity plus mechanism (ETCi / uncoupler / ATP synthase) and the inhibited ETC
complex. Mechanism is what Attene-Ramos cannot give: an uncoupler and a complex-I
inhibitor both collapse membrane potential but move OCR in opposite directions.

A drug mapping to several external records takes the most informative label
(MMP_down > cytotox > MMP_up > inactive; active > inactive), which is recorded in
`n_records` / `labels_seen` so the aggregation is auditable.

Output: one TSV per dataset with drug / group / subgroup / extra columns, written
to a PER-CELL directory, because the tables are keyed on that cell line's LINCS
perturbation names and the mappable compound set differs between cell lines.
Sharing them across cell lines would silently drop compounds.
"""
import argparse, collections, csv, os, sys

ATTENE_CLASS = {
    'active antagonist': 'MMP_down',
    'inconclusive antagonist (cytotoxic)': 'cytotox',
    'active agonist': 'MMP_up',
    'inactive': 'inactive',
}
PRIORITY = ['MMP_down', 'cytotox', 'MMP_up', 'inactive', 'incon']


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def build_attene(datatable, mapping, strict_sep, sep_fold):
    rows = [r for r in csv.DictReader(open(datatable))
            if (r.get('PUBCHEM_CID') or '').strip().isdigit()]
    by_cid = collections.defaultdict(list)
    for r in rows:
        by_cid[r['PUBCHEM_CID']].append(r)

    per_drug = collections.defaultdict(lambda: {'labels': [], 'cids': set(),
                                                'pot': [], 'eff': []})
    for m in csv.DictReader(open(mapping)):
        drug = m['pertname_dataset']
        d = per_drug[drug]
        d['cids'].add(m['val_id'])
        for r in by_cid.get(m['val_id'], []):
            lab = ATTENE_CLASS.get(r['Activity Summary'], 'incon')
            if lab == 'MMP_down' and strict_sep:
                rp, vp = fnum(r['Ratio Potency (uM)']), fnum(r['Viability Potency (uM)'])
                if rp is not None and vp is not None and vp < sep_fold * rp:
                    lab = 'cytotox'
            d['labels'].append(lab)
            if lab in ('MMP_down', 'MMP_up'):
                p, e = fnum(r['Ratio Potency (uM)']), fnum(r['Ratio Efficacy (%)'])
                if p is not None:
                    d['pot'].append(p)
                if e is not None:
                    d['eff'].append(e)

    out = []
    for drug, d in sorted(per_drug.items()):
        if not d['labels']:
            continue
        group = next((g for g in PRIORITY if g in d['labels']), 'incon')
        out.append({
            'drug': drug, 'group': group, 'subgroup': '',
            'n_records': len(d['labels']),
            'labels_seen': ';'.join(sorted(set(d['labels']))),
            'external_ids': ';'.join(sorted(d['cids'])),
            'potency_uM': '%.4g' % min(d['pot']) if d['pot'] else '',
            'efficacy_pct': '%.4g' % max(d['eff'], key=abs) if d['eff'] else '',
        })
    return out


def build_hallinger(annotated, mapping, potency):
    ann = {r['DTXSID']: r for r in csv.DictReader(open(annotated)) if r['DTXSID']}
    ac50 = {}
    if potency and os.path.exists(potency):
        for r in csv.DictReader(open(potency)):
            if r.get('hitc') == '1' and r.get('AC50_uM'):
                key = (r['dsstox_substance_id'], r['aenm'])
                ac50[key] = r['AC50_uM']

    per_drug = collections.defaultdict(lambda: {'ids': set()})
    for m in csv.DictReader(open(mapping)):
        per_drug[m['pertname_dataset']]['ids'].add(m['val_id'])

    out = []
    for drug, d in sorted(per_drug.items()):
        recs = [ann[i] for i in d['ids'] if i in ann]
        if not recs:
            continue
        group = 'active' if any(r['tier1_activity'] == 'active' for r in recs) else 'inactive'
        act = [r for r in recs if r['tier1_activity'] == 'active']
        mech = ';'.join(sorted({r['proposed_mechanism'] for r in act if r['proposed_mechanism']}))
        cplx = ';'.join(sorted({r['inhibited_ETC_complex'] for r in act
                                if r['inhibited_ETC_complex'] not in ('', 'nd')}))
        basal = [ac50[(r['DTXSID'], 'NCCT_MITO_basal_resp_rate_OCR_dn')] for r in act
                 if (r['DTXSID'], 'NCCT_MITO_basal_resp_rate_OCR_dn') in ac50]
        out.append({
            'drug': drug, 'group': group, 'subgroup': mech,
            'n_records': len(recs),
            'labels_seen': ';'.join(sorted({r['tier1_activity'] for r in recs})),
            'external_ids': ';'.join(sorted(d['ids'])),
            'inhibited_complex': cplx,
            'basal_OCR_dn_AC50_uM': basal[0] if basal else '',
            'EFA_activity': ';'.join(sorted({r['EFA_activity'] for r in act if r['EFA_activity']})),
        })
    return out


def write(rows, path, extra):
    cols = ['drug', 'group', 'subgroup', 'n_records', 'labels_seen', 'external_ids'] + extra
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        w.writerows(rows)
    counts = collections.Counter(r['group'] for r in rows)
    print('%s: %d drugs %s' % (os.path.basename(path), len(rows), dict(counts)), file=sys.stderr)


def main():
    repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--attene-datatable',
                    default=os.path.join(repo, 'data/validation/attene_ramos/AID_720637_datatable.csv'))
    ap.add_argument('--cell', default='HEPG2',
                    help='LINCS cell line whose perturbations the mapping used')
    ap.add_argument('--attene-mapping', default=None,
                    help='default: data/validation/attene_ramos/mapping_<cell>/mapped_drugs.csv')
    ap.add_argument('--hallinger-annotated',
                    default=os.path.join(repo, 'data/validation/hallinger/hallinger_tier1_all1042_annotated.csv'))
    ap.add_argument('--hallinger-mapping', default=None,
                    help='default: data/validation/hallinger/mapping_<cell>/mapped_drugs.csv')
    ap.add_argument('--hallinger-potency',
                    default=os.path.join(repo, 'data/validation/hallinger/hallinger_tier2_potency.csv'))
    ap.add_argument('--strict-sep', action='store_true',
                    help='demand MMP potency at least --sep-fold below viability potency')
    ap.add_argument('--sep-fold', type=float, default=3.0)
    ap.add_argument('--out-dir', required=True,
                    help='per-cell phenotype directory, e.g. '
                         '<scratch>/results/mc_phenotype_assoc/<cell>/phenotypes')
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # resolve the per-cell mappings, falling back to the original HEPG2 layout
    def pick(explicit, *candidates):
        if explicit:
            return explicit
        for c in candidates:
            if os.path.exists(c):
                return c
        return candidates[0]
    args.attene_mapping = pick(
        args.attene_mapping,
        os.path.join(repo, 'data/validation/attene_ramos/mapping_%s/mapped_drugs.csv' % args.cell),
        os.path.join(repo, 'data/validation/attene_ramos/mapped_drugs.csv'))
    args.hallinger_mapping = pick(
        args.hallinger_mapping,
        os.path.join(repo, 'data/validation/hallinger/mapping_%s/mapped_drugs.csv' % args.cell),
        os.path.join(repo, 'data/validation/hallinger/mapping/mapped_drugs.csv'))
    for lab, f in (('attene', args.attene_mapping), ('hallinger', args.hallinger_mapping)):
        if not os.path.exists(f):
            raise SystemExit(
                'no %s mapping for cell=%s (%s).\nRun:  cell=%s bash '
                'scripts/benchmarking/map_%s_lincs_perts.sh'
                % (lab, args.cell, f, args.cell, lab if lab != 'attene' else 'attene'))
    print('cell=%s\n  attene mapping   : %s\n  hallinger mapping: %s'
          % (args.cell, args.attene_mapping, args.hallinger_mapping), file=sys.stderr)

    a = build_attene(args.attene_datatable, args.attene_mapping, args.strict_sep, args.sep_fold)
    write(a, os.path.join(args.out_dir, 'attene_ramos_mmp.tsv'), ['potency_uM', 'efficacy_pct'])
    h = build_hallinger(args.hallinger_annotated, args.hallinger_mapping, args.hallinger_potency)
    write(h, os.path.join(args.out_dir, 'hallinger_ocr.tsv'),
          ['inhibited_complex', 'basal_OCR_dn_AC50_uM', 'EFA_activity'])
    mech = collections.Counter(r['subgroup'] for r in h if r['group'] == 'active')
    print('hallinger mechanisms: %s' % dict(mech), file=sys.stderr)


if __name__ == '__main__':
    main()
