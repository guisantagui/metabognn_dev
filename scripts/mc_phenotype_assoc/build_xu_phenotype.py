#!/usr/bin/env python3
"""Build a phenotype table from the Xu et al. AC16 cardiomyocyte MPP screen.

Source: `data/validation/xu/hp6c00052_si_002.xlsx`, sheet `Excel Table S10`,
"Experimental validation of model-predicted cardiotoxic compounds using the MPP
assay" (doi 10.1021/EHP.6c00052). 247 compounds with an assay `class`, AC50,
efficacy and Tox21 curve rank, measured in **AC16 cardiomyocytes**.

**This is a falsification test, not a cell-matched one** (report 09, 8.2 and 10, 4).
AC16 is not in LINCS at all, so no AC16 labels exist and none can be made. The
phenotype is **membrane potential - the same readout as Attene-Ramos** - measured
by a different group in a different, mitochondria-rich cell type. Scoring it
against HEPG2 labels therefore asks: do HepG2 MetChange labels predict membrane
potential in a cell line they know nothing about? A **positive** result means the
signal travels with the compounds, which would revive the chemotype explanation
and damage report 09, section 3. A **null** is weak evidence either way at this n.
Run it because it can only hurt the current conclusion.

Group definitions, chosen to be **directly comparable to the Attene-Ramos table**
so the same sets, contrasts and controls apply and the numbers can be read side by
side with report 09, section 3:

  MMP_down   class == 'active'. All 61 actives have a NEGATIVE curve rank
             (median -3.33) and NEGATIVE efficacy (median -53.6%), i.e. every
             active in this screen is a membrane-potential DECREASER - verified,
             not assumed, so the group is directionally pure.
  inactive   class == 'inactive'. All 123 have curve rank exactly 0.
  (dropped)  class == 'inconclusive', 63 compounds, mostly positive curve rank.

Note a limitation relative to the Naia table: **Xu reports no viability readout**,
so the only available confounder control is the `n_sig` responsiveness proxy, not
a measured cytotoxicity as in report 09, 9.2. Membrane-potential loss and cell
death are not independent, so an association here cannot be separated from
viability the way the Naia one could.

SMILES are canonicalised to an InChIKey with RDKit and matched to
`lincs_unique_perts_<CELL>.tsv` on the first block, which ignores salt and stereo
and is therefore a floor; `map_smiles.R` would recover more.
"""
import argparse, collections, csv, os, re, sys

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def norm(s):
    s = (s or '').strip().lower()
    s = (s.replace('’', "'").replace('′', "'")
          .replace('–', '-').replace('−', '-'))
    s = re.sub(r'[\s,;]+', '', s)
    return re.sub(r"[\(\)\[\]\{\}'\"`.]", '', s)


def num(v):
    """'-3.33+-0.58' or '  2.65+- 0.30' -> the mean as a float."""
    if isinstance(v, (int, float)):
        return float(v)
    if v is None:
        return None
    m = re.match(r'^\s*(-?[0-9.]+)', str(v).replace('−', '-'))
    return float(m.group(1)) if m else None


def read_xu(path, sheet):
    from openpyxl import load_workbook
    wb = load_workbook(path, read_only=True, data_only=True)
    if sheet not in wb.sheetnames:
        raise SystemExit('sheet %r not in %s (have: %s)' % (sheet, path, wb.sheetnames))
    rows = [r for r in wb[sheet].iter_rows(values_only=True)]
    wb.close()
    hi = next((i for i, r in enumerate(rows)
               if r and 'SMILES' in [str(x) for x in r]), None)
    if hi is None:
        raise SystemExit('no header row containing SMILES in sheet %r' % sheet)
    hdr = [str(x) if x is not None else '' for x in rows[hi]]
    out = []
    for r in rows[hi + 1:]:
        if not r or not any(r):
            continue
        d = dict(zip(hdr, r))
        if not d.get('SMILES'):
            continue
        out.append({'name': str(d.get('Sample.Name') or '').strip(),
                    'smiles': str(d['SMILES']).strip(),
                    'cls': str(d.get('class') or '').strip().lower(),
                    'tox21': str(d.get('TOX21_ID') or '').strip(),
                    'ac50': num(d.get('AC50')),
                    'eff': num(d.get('Efficacy')),
                    'rank': num(d.get('Curve.rank'))})
    return out


def lincs_index(path):
    by_ik, by_nm = collections.defaultdict(list), {}
    with open(path, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            nm = (r.get('pertname') or '').strip()
            if not nm:
                continue
            ik = (r.get('inchi_key') or '').strip()
            if ik and ik != 'NA':
                by_ik[ik.split('-')[0]].append(nm)
            by_nm.setdefault(norm(nm), nm)
    return by_ik, by_nm


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--xu', default=os.path.join(REPO, 'data/validation/xu/hp6c00052_si_002.xlsx'))
    ap.add_argument('--sheet', default='Excel Table S10')
    ap.add_argument('--cell', default='HEPG2', help='which LINCS panel to map drug names against')
    ap.add_argument('--lincs-perts', default='')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    perts = args.lincs_perts or os.path.join(REPO, 'data/LINCS/lincs_unique_perts_%s.tsv' % args.cell)
    if not os.path.exists(perts):
        raise SystemExit('no LINCS panel file %s' % perts)

    from rdkit import Chem, RDLogger
    from rdkit.Chem.inchi import MolToInchiKey
    RDLogger.DisableLog('rdApp.*')

    rows = read_xu(args.xu, args.sheet)
    by_ik, by_nm = lincs_index(perts)
    stats = collections.Counter()

    # verify the directional purity the docstring asserts, rather than trusting it
    act = [r for r in rows if r['cls'] == 'active']
    bad = [r for r in act if (r['rank'] or 0) >= 0]
    print('class counts: %s' % dict(collections.Counter(r['cls'] for r in rows)), file=sys.stderr)
    print('actives with a non-negative curve rank: %d of %d%s'
          % (len(bad), len(act), ' - group is NOT directionally pure' if bad else ' (all decreasers)'),
          file=sys.stderr)

    best = {}
    for r in rows:
        if r['cls'] not in ('active', 'inactive'):
            stats['dropped_' + (r['cls'] or 'blank')] += 1
            continue
        mol = Chem.MolFromSmiles(r['smiles'])
        ik = MolToInchiKey(mol) if mol is not None else ''
        if not ik:
            stats['smiles_unparsed'] += 1
        hits, via = (by_ik.get(ik.split('-')[0], []), 'inchikey') if ik else ([], '')
        if not hits:
            nm = by_nm.get(norm(r['name']))
            hits, via = ([nm], 'name') if nm else ([], '')
        if not hits:
            stats['no_lincs'] += 1
            continue
        stats['matched_' + via] += 1
        group = 'MMP_down' if r['cls'] == 'active' else 'inactive'
        for pertname in hits:
            prev = best.get(pertname)
            # a LINCS drug matching both an active and an inactive row is a conflict;
            # keep the active call and count it, rather than letting order decide
            if prev is None:
                best[pertname] = dict(r, group=group, via=via, n_records=1)
            elif prev['group'] != group:
                stats['conflicting_calls'] += 1
                if group == 'MMP_down':
                    best[pertname] = dict(r, group=group, via=via,
                                          n_records=prev['n_records'] + 1)
                else:
                    prev['n_records'] += 1
            else:
                prev['n_records'] += 1

    cols = ['drug', 'group', 'subgroup', 'n_records', 'labels_seen', 'external_ids',
            'potency_uM', 'efficacy_pct', 'curve_rank', 'match_via']
    n = collections.Counter()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or '.', exist_ok=True)
    with open(args.out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(cols)
        for pertname, r in sorted(best.items()):
            n[r['group']] += 1
            f = lambda v: '' if v is None else '%.4g' % v
            w.writerow([pertname, r['group'], '', r['n_records'], r['cls'], r['tox21'],
                        f(r['ac50']), f(r['eff']), f(r['rank']), r['via']])
    print('mapping to the %s panel: %s' % (args.cell, dict(stats)), file=sys.stderr)
    print('groups written: %s' % dict(n), file=sys.stderr)
    print('wrote %s' % args.out, file=sys.stderr)


if __name__ == '__main__':
    main()
