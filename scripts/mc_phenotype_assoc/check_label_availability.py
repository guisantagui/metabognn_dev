#!/usr/bin/env python3
"""Which cell lines actually have the limma label matrices this analysis needs?

Checked at runtime rather than hardcoded, because the answer differs between the
laptop and the cluster and changes as label runs finish.

Per cell line and per arm:
  REQUIRED   log2FC.csv   - the `value` and `sign` modalities read it
  OPTIONAL   qvals.csv    - needed by `--mode sig`, and used for n_sig (hence the
                            responsiveness stratification) when fit_info_nsig.csv
                            is absent
  OPTIONAL   fit_info_nsig.csv - gives n_sig and the parsed drug/dose/time
                            directly; without it both are derived (drug/dose/time
                            from the column names, n_sig from the q-values)

A cell line counts as runnable for an arm when log2FC.csv is present and
non-empty. `--list {mc,expr,both}` prints just the runnable cell names, one per
line, for a submit script to consume.
"""
import argparse, os, sys


def arm_paths(base, cell, arm):
    sub = 'metchange/rank' if arm == 'mc' else 'expr'
    return os.path.join(base, 'results', 'labels', cell, sub)


def inspect(d):
    out = {}
    for f in ('log2FC.csv', 'qvals.csv', 'fit_info_nsig.csv'):
        p = os.path.join(d, f)
        out[f] = os.path.isfile(p) and os.path.getsize(p) > 0
    out['_dir'] = os.path.isdir(d)
    return out


def main():
    repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--labels-base', default=os.environ.get('labels_base', repo),
                    help='project root holding results/labels/<cell>/... '
                         '(scratch on the cluster)')
    ap.add_argument('--cell-list', default=os.path.join(repo, 'data/LINCS/core_cell_lines.tsv'))
    ap.add_argument('--cells', default='', help='space/comma separated, overrides --cell-list')
    ap.add_argument('--list', choices=['mc', 'expr', 'both'],
                    help='print only the runnable cell names (both = runnable for either arm)')
    args = ap.parse_args()

    if args.cells.strip():
        cells = [c for c in args.cells.replace(',', ' ').split() if c]
    else:
        with open(args.cell_list) as fh:
            cells = [l.strip().strip('"') for l in fh if l.strip()
                     and not l.startswith('#')]
        # tolerate a header line
        if cells and cells[0].lower() in ('cell', 'cell_line', 'cell_id'):
            cells = cells[1:]

    rows = []
    for c in cells:
        r = {'cell': c}
        for arm in ('mc', 'expr'):
            d = arm_paths(args.labels_base, c, arm)
            i = inspect(d)
            r[arm] = i
            r[arm + '_ok'] = i['log2FC.csv']
        rows.append(r)

    if args.list:
        for r in rows:
            ok = (r['mc_ok'] or r['expr_ok']) if args.list == 'both' else r[args.list + '_ok']
            if ok:
                print(r['cell'])
        return

    print('labels base: %s\n' % args.labels_base)
    hdr = '%-9s %-28s %-28s %s' % ('cell', 'metchange/rank', 'expr', 'runnable')
    print(hdr); print('-' * len(hdr))
    def cell_str(i):
        if not i['_dir']:
            return 'no directory'
        bits = []
        for f, tag in (('log2FC.csv', 'log2FC'), ('qvals.csv', 'qvals'),
                       ('fit_info_nsig.csv', 'nsig')):
            bits.append(tag if i[f] else '-' * len(tag))
        return ' '.join(bits)
    n_mc = n_ex = 0
    for r in rows:
        which = []
        if r['mc_ok']:
            which.append('mc'); n_mc += 1
        if r['expr_ok']:
            which.append('expr'); n_ex += 1
        print('%-9s %-28s %-28s %s' % (r['cell'], cell_str(r['mc']), cell_str(r['expr']),
                                       ','.join(which) if which else 'NO'))
    print('\nrunnable: %d cell line(s) for the MetChange arm, %d for the expression arm'
          % (n_mc, n_ex))
    missing = [r['cell'] for r in rows if not (r['mc_ok'] or r['expr_ok'])]
    if missing:
        print('no usable labels: %s' % ', '.join(missing))
    nonsig = [r['cell'] for r in rows if r['mc_ok']
              and not (r['mc']['qvals.csv'] or r['mc']['fit_info_nsig.csv'])]
    if nonsig:
        print('WARNING no q-values and no fit_info_nsig for the MetChange arm of: %s'
              % ', '.join(nonsig))
        print('  -> --mode sig unavailable and n_sig unknown, so the responsiveness')
        print('     stratification (report 08, 3) cannot run for those cell lines.')


if __name__ == '__main__':
    main()
