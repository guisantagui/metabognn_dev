#!/usr/bin/env python3
"""Parse the Hallinger et al. 2020 supplemental tables and map them to LINCS.

The Dryad deposit names only the 270 tier-2 samples; the ~780 remaining tier-1
spids are EPA source-plate wells (EPAPLT0033/0034/0035) with no public plate map.
The plate map is not needed: Supplemental Table 2 of the paper lists all 1,042
tier-1 chemicals by name with their activity calls, and a second supplemental
table characterises the 243 actives (per-aeid calls, proposed mechanism, electron
flow assay, inhibited ETC complex).

Both tables are PDFs, hosted on EPA Clowder and linked from the EPA figshare
record doi:10.23645/epacomptox.12123735:
  SC summary (19 pp, all 1,042)   file id 66969672e4b0a7c65d1a38a1
  Study summary (26 pp, 243 act.) file id 668451dae4b0a7c65d175744
  QC summary                      file id 66969672e4b0a7c65d1a3899
Download with: curl -L -o <name>.pdf \
  https://clowder.edap-cluster.com/api/files/<file id>/blob
(HEAD 404s on that endpoint, GET works.)

Outputs, in data/validation/hallinger/:
  hallinger_tier1_suppT2.csv            1,042 x (name, conc, activity)
  hallinger_tier2_mechanism.csv         243 x (per-aeid calls, mechanism, complex)
  hallinger_tier1_all1042_annotated.csv the above + DSSTox structures + LINCS HepG2 match

Needs pypdf. Name -> structure uses the dated DSSTox CCD dump (PREFERRED_NAME and
IDENTIFIER, i.e. synonyms); LINCS matching here is on the first InChIKey block,
which is a floor - run map_smiles.R on the SMILES column for the proper mapping
(salt stripping, stereo handling).
"""
import argparse, collections, csv, glob, os, re, sys

csv.field_size_limit(10 ** 9)

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
HALL = os.path.join(REPO, 'data', 'validation', 'hallinger')

MECH_COLS = ['proposed_mechanism', 'EFA_activity', 'inhibited_ETC_complex',
             'decr_basal_aeid2442', 'decr_maximal_aeid2444',
             'incr_basal_aeid2443', 'incr_inhibited_aeid2447']


def demojibake(s):
    """Repair UTF-8-read-as-latin1 artefacts in the PDF text (e.g. "Nâ€²" for "N'")."""
    if 'Ã' in s or 'â' in s:
        try:
            s = s.encode('latin1').decode('utf-8')
        except (UnicodeEncodeError, UnicodeDecodeError):
            pass
    return s


def norm(s):
    s = demojibake(s)
    s = s.strip().lower().replace('’', "'").replace('‘', "'")
    s = s.replace('′', "'").replace('″', "'")
    s = s.replace('–', '-').replace('−', '-')
    s = re.sub(r'[\s,;]+', '', s)
    return re.sub(r"[\(\)\[\]\{\}'\"`.]", '', s)


def pdf_lines(path):
    import pypdf
    out = []
    for page in pypdf.PdfReader(path).pages:
        out += [l.strip() for l in page.extract_text().split('\n') if l.strip()]
    return out


def parse_tier1(path):
    """Supplemental Table 2, SC summary: chemical, test concentration, activity."""
    pat = re.compile(r'^(.*?)\s+([0-9]+(?:\.[0-9]+)?)\s+(active|inactive)$', re.I)
    rows = []
    for line in pdf_lines(path):
        if line.lower().startswith('chemical '):
            continue
        m = pat.match(line)
        if m:
            rows.append([demojibake(m.group(1).strip()), m.group(2), m.group(3).lower()])
    return rows


def parse_mechanism(path):
    """Study summary: the 243 tier-1 actives with their tier-2 / EFA characterisation."""
    pat = re.compile(r'^(.+?)\s+([0-9]+(?:\.[0-9]+)?)\s+(active)\s+'
                     r'(active|inactive)\s+(active|inactive)\s+(active|inactive)\s+(active|inactive)\s+'
                     r'(ETCi|Uncoupler|ATP Synthase|Redox|Inactive|Interference|[A-Za-z /]+?)'
                     r'(?:\s+(active|inactive)\s+(.+))?$')
    rows = []
    for line in pdf_lines(path):
        m = pat.match(line)
        if not m:
            continue
        g = list(m.groups())
        # PDF column order is 2447, 2443, 2444, 2442; emit it as 2442, 2444, 2443, 2447
        rows.append([demojibake(g[0].strip()), g[1], g[2], g[6], g[5], g[4], g[3],
                     g[7].strip(), g[8] or '', g[9] or ''])
    return rows


def dsstox_lookup(names, dump_dir):
    """normalised name -> (DTXSID, PREFERRED_NAME, INCHIKEY, SMILES, via)."""
    want = set(names)
    hit = {}
    files = sorted(glob.glob(os.path.join(dump_dir, 'DSSToxCCDdump[0-9]*.csv')))
    if not files:
        print('no DSSTox dump shards in %s - structures left empty' % dump_dir, file=sys.stderr)
        return hit
    rank = {'preferred': 0, 'synonym': 1}
    for fp in files:
        with open(fp, newline='', encoding='utf-8', errors='replace') as fh:
            for row in csv.DictReader(fh):
                # IDENTIFIER is a pipe-separated synonym list (CASRNs and names), so split it
                cands = [(row.get('PREFERRED_NAME') or '', 'preferred')]
                cands += [(p, 'synonym') for p in (row.get('IDENTIFIER') or '').split('|')]
                for val, via in cands:
                    if not val.strip():
                        continue
                    key = norm(val)
                    if key not in want:
                        continue
                    prev = hit.get(key)
                    if prev is None or rank[via] < rank[prev[4]]:
                        hit[key] = (row['DTXSID'], row['PREFERRED_NAME'],
                                    row['INCHIKEY'], row['SMILES'], via)
    return hit


def lincs_index(path):
    idx = collections.defaultdict(list)
    with open(path, newline='') as fh:
        for row in csv.DictReader(fh, delimiter='\t'):
            key = (row.get('inchi_key') or '').strip()
            if key and key != 'NA':
                idx[key.split('-')[0]].append(row)
    return idx


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--sc-pdf', default=os.path.join(HALL, 'supplemental', 'SuppTable2_SC.summ.table.pdf'))
    ap.add_argument('--study-pdf', default=os.path.join(HALL, 'supplemental', 'SuppTable2_Study_summary.pdf'))
    ap.add_argument('--dsstox-dir', default=os.path.join(HALL, 'DSSTox_CCD_dump_12092025'))
    ap.add_argument('--lincs', default=os.path.join(REPO, 'data', 'LINCS', 'lincs_unique_perts_HEPG2.tsv'))
    ap.add_argument('--out-dir', default=HALL)
    args = ap.parse_args()

    tier1 = parse_tier1(args.sc_pdf)
    print('tier 1: %d chemicals (%s)' % (
        len(tier1), dict(collections.Counter(r[2] for r in tier1))))
    with open(os.path.join(args.out_dir, 'hallinger_tier1_suppT2.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['chnm', 'test_conc_uM', 'activity'])
        w.writerows(tier1)

    mech_rows = parse_mechanism(args.study_pdf)
    print('tier 1 actives characterised: %d (%s)' % (
        len(mech_rows), dict(collections.Counter(r[7] for r in mech_rows))))
    mech_header = ['chnm', 'test_conc_uM', 'tier1_activity', 'decr_basal_aeid2442',
                   'decr_maximal_aeid2444', 'incr_basal_aeid2443', 'incr_inhibited_aeid2447',
                   'proposed_mechanism', 'EFA_activity', 'inhibited_ETC_complex']
    with open(os.path.join(args.out_dir, 'hallinger_tier2_mechanism.csv'), 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(mech_header)
        w.writerows(mech_rows)
    mech = {norm(r[0]): dict(zip(mech_header, r)) for r in mech_rows}

    struct = dsstox_lookup([norm(r[0]) for r in tier1], args.dsstox_dir)
    print('names resolved to a structure: %d of %d' % (len(struct), len(tier1)))

    lincs = lincs_index(args.lincs) if os.path.exists(args.lincs) else {}
    matched = collections.Counter()
    out_path = os.path.join(args.out_dir, 'hallinger_tier1_all1042_annotated.csv')
    with open(out_path, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['chnm', 'test_conc_uM', 'tier1_activity', 'DTXSID', 'PREFERRED_NAME',
                    'INCHIKEY', 'SMILES', 'name_match_via'] + MECH_COLS
                   + ['lincs_hepg2_pert_id', 'lincs_hepg2_pertname', 'n_lincs_matches'])
        for name, conc, activity in tier1:
            s = struct.get(norm(name), ('', '', '', '', ''))
            m = mech.get(norm(name), {})
            hits = lincs.get(s[2].split('-')[0], []) if s[2] else []
            if hits:
                matched[activity] += 1
            w.writerow([name, conc, activity, s[0], s[1], s[2], s[3], s[4]]
                       + [m.get(c, '') for c in MECH_COLS]
                       + ['|'.join(h['pert_id'] for h in hits),
                          '|'.join(h['pertname'] for h in hits), len(hits)])
    print('LINCS HepG2 matched (InChIKey first block, floor): %d (%s)'
          % (sum(matched.values()), dict(matched)))
    print('wrote', out_path)


if __name__ == '__main__':
    main()
