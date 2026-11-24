#!/usr/bin/env python3
"""Build the PRISM phenotype table for the within-drug, between-cell differential.

**Run this on the laptop, where the PRISM download lives, and upload the single
small TSV it writes.** The PRISM files are ~800 MB in total but the table is a few
hundred kB, so nothing large has to move in either direction: the label matrices
stay on the cluster and this table goes to them.

Source: `data/validation/prism/secondary-screen-dose-response-curve-parameters.csv`
(PRISM Repurposing secondary screen, depmap.org/repurposing). 1,449 BRD compounds x
482 cell lines, 8-step 4-fold dilution from 10 uM. `broad_id` is a BRD id and
`ccle_name` names the line, so mapping needs no structure work - only
`pert_id -> pertname` from `data/LINCS/lincs_unique_perts.tsv`, which is global
rather than per cell line.

WHY THE TABLE LOOKS LIKE THIS (report 10, 1.1-1.2). The PRISM phenotype is
viability, and viability is the confounder that dominated reports 07-09 - report 09,
9.6 turns on it, with the expression association dying under a cytotoxicity match
while MetChange's survived. **You cannot match on cytotoxicity when cytotoxicity is
the phenotype**, so the per-cell "do sensitive drugs score differently from
resistant ones" test is close to uninterpretable: any generally bioactive compound
lowers AUC and also produces large MetChange signatures in every set.

The design that works is the **differential**: per drug, correlate
`AUC(X) - AUC(Y)` against `setscore(X labels) - setscore(Y labels)`. Differencing
removes the drug's general potency by construction, so a compound that is simply
toxic everywhere contributes zero to both sides and only the *selectivity* question
remains. `prism_differential.py` computes it from the per-cell scores.

This script therefore writes a table whose **groups exist only to satisfy
`mc_phenotype_assoc.py`'s interface**, which requires a case and a control group.
They are a cell-INDEPENDENT median split on each drug's mean AUC across the covered
cell lines, so that (a) every drug is scored in every cell line, which the
differential needs, and (b) the split is identical in all runs. The per-cell
contrast those groups produce is the confounded naive test, and
`prism_differential.py` prints it only as a foil, never as the result.

The per-cell AUCs are carried as `auc_<CELL>` columns; they are the phenotype side
of the differential.

PRE-REGISTERED DIRECTION, written down before the test is run (report 10, 1.2):
lower predicted production capacity for biomass precursors in X than in Y should go
with MORE growth inhibition in X, i.e. LOWER AUC in X. Both differences are then
negative together, so the expected correlation is **positive**.
"""
import argparse, collections, csv, os, sys

csv.field_size_limit(10 ** 9)
REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# LINCS core cell lines present in the PRISM secondary screen, verified 2026-11-24.
# HEPG2, VCAP and HCC515 are ABSENT (PRISM has 14 liver lines, none of them HEPG2);
# HA1E is absent by construction, being an immortalised normal line rather than a
# CCLE cancer line. See report 10, section 1.
CCLE = {'A375': 'A375_SKIN', 'A549': 'A549_LUNG', 'HT29': 'HT29_LARGE_INTESTINE',
        'MCF7': 'MCF7_BREAST', 'PC3': 'PC3_PROSTATE'}


def brd(s):
    """'BRD-K12345678-001-02-9' -> 'BRD-K12345678', the LINCS pert_id form."""
    return '-'.join(str(s).split('-')[:2])


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--prism', default=os.path.join(
        REPO, 'data/validation/prism/secondary-screen-dose-response-curve-parameters.csv'))
    ap.add_argument('--lincs-perts', default=os.path.join(REPO, 'data/LINCS/lincs_unique_perts.tsv'),
                    help='global pert_id -> pertname map; NOT a per-cell panel file')
    ap.add_argument('--cells', default=' '.join(sorted(CCLE)),
                    help='LINCS core names to include (must be PRISM-covered)')
    ap.add_argument('--metric', default='auc', choices=['auc', 'ic50', 'ec50'],
                    help="phenotype column; 'auc' is the one to use - it is defined for "
                         'every drug x line, while ic50/ec50 are missing for inactives')
    ap.add_argument('--require-all', action='store_true', default=True,
                    help='keep only drugs measured in every requested cell line (default); '
                         'the common-drug restriction of report 09, section 2')
    ap.add_argument('--allow-partial', dest='require_all', action='store_false',
                    help='keep drugs missing in some lines (the differential then uses '
                         'different drug sets per cell pair - read report 09, section 2 first)')
    ap.add_argument('--min-r2', type=float, default=0.0,
                    help='drop curve fits below this r2 (0 = keep all)')
    ap.add_argument('--out', required=True)
    args = ap.parse_args()

    cells = [c for c in args.cells.replace(',', ' ').split() if c]
    unknown = [c for c in cells if c not in CCLE]
    if unknown:
        raise SystemExit('not PRISM-covered (see report 10, section 1): %s\ncovered: %s'
                         % (', '.join(unknown), ', '.join(sorted(CCLE))))
    want = {CCLE[c]: c for c in cells}

    # pert_id -> pertname, the `drug` field of a label dir's fit_info_nsig.csv
    id2name = {}
    with open(args.lincs_perts, newline='') as fh:
        for r in csv.DictReader(fh, delimiter='\t'):
            pid, nm = (r.get('pert_id') or '').strip(), (r.get('pertname') or '').strip()
            if pid and nm:
                id2name.setdefault(pid, nm)
    print('LINCS pert_id -> pertname: %d entries' % len(id2name), file=sys.stderr)

    vals = collections.defaultdict(dict)        # pertname -> cell -> [metric,...]
    moa, stats = {}, collections.Counter()
    with open(args.prism, newline='') as fh:
        for r in csv.DictReader(fh):
            cell = want.get(r['ccle_name'])
            if cell is None:
                continue
            stats['rows_in_wanted_cells'] += 1
            name = id2name.get(brd(r['broad_id']))
            if name is None:
                stats['no_lincs_pert_id'] += 1
                continue
            try:
                v = float(r[args.metric])
            except (ValueError, KeyError):
                stats['no_%s' % args.metric] += 1
                continue
            if args.min_r2:
                try:
                    if float(r['r2']) < args.min_r2:
                        stats['low_r2'] += 1
                        continue
                except (ValueError, KeyError):
                    pass
            vals[name].setdefault(cell, []).append(v)
            if r.get('moa'):
                moa.setdefault(name, r['moa'])
    print('PRISM: %s' % dict(stats), file=sys.stderr)

    # several PRISM broad_ids can collapse to one LINCS pertname; average them
    per = {n: {c: sum(v) / len(v) for c, v in d.items()} for n, d in vals.items()}
    if args.require_all:
        per = {n: d for n, d in per.items() if len(d) == len(cells)}
    if not per:
        raise SystemExit('no drugs left after mapping and filtering')

    # cell-INDEPENDENT median split, purely to satisfy the case/control interface
    mean_auc = {n: sum(d.values()) / len(d) for n, d in per.items()}
    order = sorted(mean_auc, key=lambda n: mean_auc[n])
    mid = len(order) // 2
    group = {n: ('sensitive' if i < mid else 'resistant') for i, n in enumerate(order)}

    cols = (['drug', 'group', 'subgroup', 'n_records', 'labels_seen', 'external_ids',
             'mean_%s' % args.metric] + ['%s_%s' % (args.metric, c) for c in cells] + ['moa'])
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or '.', exist_ok=True)
    with open(args.out, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(cols)
        for n in order:
            d = per[n]
            w.writerow([n, group[n], '', len(vals[n]), args.metric, '',
                        '%.5g' % mean_auc[n]]
                       + ['' if c not in d else '%.5g' % d[c] for c in cells]
                       + [moa.get(n, '')])
    import statistics as st
    spread = [st.pstdev(list(per[n].values())) for n in order if len(per[n]) > 1]
    print('\ndrugs written: %d  (cells: %s)' % (len(order), ', '.join(cells)), file=sys.stderr)
    print('groups (interface only, cell-independent median split on mean %s): %s'
          % (args.metric, dict(collections.Counter(group.values()))), file=sys.stderr)
    if spread:
        spread.sort()
        print('between-cell %s spread per drug: median %.3f, 75th %.3f, max %.3f'
              % (args.metric, spread[len(spread) // 2], spread[int(.75 * len(spread))], spread[-1]),
              file=sys.stderr)
        print('drugs with spread > 0.1 (the informative ones for the differential): %d'
              % sum(1 for x in spread if x > 0.1), file=sys.stderr)
    print('wrote %s' % args.out, file=sys.stderr)
    print('\nNext: copy this ONE file to the cluster, then run the launcher there:', file=sys.stderr)
    print('  scp %s <host>:<scratch_proj_dir>/results/mc_phenotype_assoc/phenotypes_ALL/'
          % args.out, file=sys.stderr)
    print('  cd scripts/mc_phenotype_assoc && sbatch launch_prism_slurm.sh', file=sys.stderr)


if __name__ == '__main__':
    main()
