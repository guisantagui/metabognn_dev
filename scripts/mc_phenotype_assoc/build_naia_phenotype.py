#!/usr/bin/env python3
"""Build a phenotype table from the Naia et al. SH-SY5Y ATP screen.

Naia et al., BMC Biology 2021 (PMC7989211) screened ~1,200 FDA-approved compounds
in differentiated SH-SY5Y cells for cellular ATP, reporting a z-score for ATP and
one for cytotoxicity per compound. Source file:
`data/validation/naia/12915_2021_979_MOESM2_ESM.xlsx`.

**This is NOT the cell-matched test.** LINCS has only 224 perturbations in SHSY5Y
and just 10 of Naia's compounds are among them (report 09, 8.1), so SH-SY5Y labels
are impossible. What this table is for is the converse design: score the phenotype
against **another** cell line's labels, usually HEPG2, and ask whether labels from
one cell line predict a phenotype measured in a different one. Read the direction
carefully - report 09, section 3 fixes the phenotype and varies the label cell
line, while this fixes the labels and varies the phenotype cell line. A positive
result here is evidence AGAINST the cell-specificity reading of report 09.

Groups, from the ATP z-score, with a symmetric dead band so the control group is
"measured and flat" rather than "everything else":
  ATP_down   z <= -(--z-cut)      ATP lost
  ATP_up     z >= +(--z-cut)      ATP gained (the paper's own finding: enhancers)
  inactive   |z| < --z-flat       measured, no ATP effect
  excluded   everything in between

Cytotoxicity is deliberately NOT used to filter. It is carried through as a column
so the downstream test can match on it, which is the better control: it is a
direct per-compound viability measurement rather than the `n_sig` responsiveness
proxy used in reports 07-09, and in this screen the two readouts are nearly
independent (Spearman(z ATP, z cytotox) = -0.083 over 1,199 compounds).

Mapping to LINCS drug names goes name/CASRN -> DSSTox -> InChIKey first block ->
`lincs_unique_perts_<CELL>.tsv` pertname, which is the `drug` field of the label
directory's `fit_info_nsig.csv`. The first block ignores salt and stereo, so this
is a floor; `map_smiles.R` would recover more.
"""
import argparse, collections, csv, glob, os, re, sys

csv.field_size_limit(10 ** 9)
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def norm(s):
    s = (s or '').strip().lower()
    s = (s.replace('’', "'").replace('‘', "'").replace('′', "'")
          .replace('″', "'").replace('–', '-').replace('−', '-'))
    s = re.sub(r'[\s,;]+', '', s)
    return re.sub(r"[\(\)\[\]\{\}'\"`.]", '', s)


def read_naia(path):
    from openpyxl import load_workbook
    wb = load_workbook(path, read_only=True, data_only=True)
    ws = wb.worksheets[0]
    it = ws.iter_rows(values_only=True)
    hdr = next(it)
    hi = {h: i for i, h in enumerate(hdr) if h}
    for need in ['Chemical name', 'CAS number', 'z-score ATP', 'z-score cytotoxicity']:
        if need not in hi:
            raise SystemExit('column %r missing from %s' % (need, path))
    out = []
    for row in it:
        name = row[hi['Chemical name']]
        if not name:
            continue

        def num(col):
            v = row[hi[col]]
            try:
                return float(v)
            except (TypeError, ValueError):
                return None
        out.append({'name': str(name).strip(),
                    'cas': str(row[hi['CAS number']] or '').strip(),
                    'moa': str(row[hi['Mechanism of action']] or '').strip() if 'Mechanism of action' in hi else '',
                    'atp': num('z-score ATP'), 'cyt': num('z-score cytotoxicity')})
    wb.close()
    return out


def dsstox_keys(rows, dump_dir):
    """(name|cas) -> InChIKey, from the dated DSSTox CCD dump."""
    want_cas = {r['cas'] for r in rows if r['cas']}
    want_nm = {norm(r['name']) for r in rows}
    files = sorted(glob.glob(os.path.join(dump_dir, 'DSSToxCCDdump*.csv')))
    if not files:
        raise SystemExit('no DSSTox dump shards in %s' % dump_dir)
    cas2ik, nm2ik = {}, {}
    for fp in files:
        with open(fp, newline='', encoding='utf-8', errors='replace') as fh:
            for row in csv.DictReader(fh):
                ik = (row.get('INCHIKEY') or '').strip()
                if not ik:
                    continue
                c = (row.get('CASRN') or '').strip()
                if c and c in want_cas:
                    cas2ik.setdefault(c, ik)
                # IDENTIFIER is a pipe-separated synonym list
                for val in [row.get('PREFERRED_NAME') or ''] + (row.get('IDENTIFIER') or '').split('|'):
                    k = norm(val)
                    if k and k in want_nm:
                        nm2ik.setdefault(k, ik)
    return cas2ik, nm2ik


def lincs_index(path):
    """InChIKey first block -> [pertname]; plus normalised pertname -> [pertname]."""
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
    ap.add_argument('--naia', default=os.path.join(REPO, 'data/validation/naia/12915_2021_979_MOESM2_ESM.xlsx'))
    ap.add_argument('--dsstox-dir', default=os.path.join(REPO, 'data/validation/hallinger/DSSTox_CCD_dump_12092025'))
    ap.add_argument('--cell', default='HEPG2', help='which LINCS panel to map drug names against')
    ap.add_argument('--lincs-perts', default='', help='default data/LINCS/lincs_unique_perts_<CELL>.tsv')
    ap.add_argument('--z-cut', type=float, default=2.0, help='|z ATP| at or beyond which a compound is a case')
    ap.add_argument('--z-flat', type=float, default=0.5, help='|z ATP| below which a compound is a control')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    perts = args.lincs_perts or os.path.join(REPO, 'data/LINCS/lincs_unique_perts_%s.tsv' % args.cell)
    if not os.path.exists(perts):
        raise SystemExit('no LINCS panel file %s - download it for %s first' % (perts, args.cell))
    rows = read_naia(args.naia)
    cas2ik, nm2ik = dsstox_keys(rows, args.dsstox_dir)
    by_ik, by_nm = lincs_index(perts)

    def group(z):
        if z is None:
            return 'excluded'
        if z <= -args.z_cut:
            return 'ATP_down'
        if z >= args.z_cut:
            return 'ATP_up'
        if abs(z) < args.z_flat:
            return 'inactive'
        return 'excluded'

    best, stats = {}, collections.Counter()
    for r in rows:
        ik = cas2ik.get(r['cas']) or nm2ik.get(norm(r['name']), '')
        hits = by_ik.get(ik.split('-')[0], []) if ik else []
        via = 'inchikey'
        if not hits:
            nm = by_nm.get(norm(r['name']))
            hits, via = ([nm], 'name') if nm else ([], '')
        stats['resolved' if ik else 'unresolved'] += 1
        if not hits:
            stats['no_lincs_%s' % args.cell] += 1
            continue
        stats['matched_via_%s' % via] += 1
        g = group(r['atp'])
        for pertname in hits:
            # one drug can match several naia rows; keep the most extreme ATP effect
            prev = best.get(pertname)
            if prev is None or abs(r['atp'] or 0) > abs(prev['atp'] or 0):
                best[pertname] = dict(r, group=g, via=via, n_records=1)
            else:
                prev['n_records'] += 1

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or '.', exist_ok=True)
    cols = ['drug', 'group', 'subgroup', 'n_records', 'labels_seen', 'external_ids',
            'z_atp', 'z_cytotox', 'moa', 'match_via']
    n = collections.Counter()
    with open(args.out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(cols)
        for pertname, r in sorted(best.items()):
            if r['group'] == 'excluded':
                n['excluded'] += 1
                continue
            n[r['group']] += 1
            w.writerow([pertname, r['group'], '', r['n_records'], r['group'], r['cas'],
                        '' if r['atp'] is None else '%.4g' % r['atp'],
                        '' if r['cyt'] is None else '%.4g' % r['cyt'],
                        r['moa'], r['via']])
    print('naia: %d rows, %d resolved to a structure' % (len(rows), stats['resolved']), file=sys.stderr)
    print('mapped to the %s panel: %d by InChIKey, %d by name, %d unmatched'
          % (args.cell, stats['matched_via_inchikey'], stats['matched_via_name'],
             stats['no_lincs_%s' % args.cell]), file=sys.stderr)
    print('groups written (z_cut=%g, z_flat=%g): %s  [%d excluded as intermediate]'
          % (args.z_cut, args.z_flat,
             {k: v for k, v in n.items() if k != 'excluded'}, n['excluded']), file=sys.stderr)
    print('wrote %s' % args.out, file=sys.stderr)


if __name__ == '__main__':
    main()
