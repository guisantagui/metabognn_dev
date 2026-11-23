#!/usr/bin/env python3
"""Cross-tabulate the group calls of two phenotype tables over their shared drugs.

Asks whether two phenotype datasets are measuring the same thing before their
association results are compared. Report 09, 9.5 uses it on the Attene-Ramos
membrane-potential table against the Naia ATP table: they share 238 drugs, but
only 9 of the 23 ATP-down drugs with an MMP call are MMP-down and 9 are
MMP-inactive, so the two readouts are correlated (an ATP-down drug is about 3.8x
more likely to be MMP-down than an inactive drug is) and far from identical. That
matters, because it is why the two associations implicating different metabolite
sets is unsurprising, and why one shared factor - chemical or stress-related -
explains both poorly.
"""
import argparse, collections, csv

try:
    from scipy import stats as sps
except Exception:                                   # noqa: BLE001
    sps = None


def read(path):
    return {r['drug']: r['group'] for r in csv.DictReader(open(path), delimiter='\t')}


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--a', required=True, help='first phenotype table')
    ap.add_argument('--b', required=True, help='second phenotype table')
    ap.add_argument('--a-name', default='A')
    ap.add_argument('--b-name', default='B')
    ap.add_argument('--a-case', default='', help='optional: case group of A, for the 2x2 odds ratio')
    ap.add_argument('--b-case', default='', help='optional: case group of B')
    ap.add_argument('--control', default='inactive')
    args = ap.parse_args()

    A, B = read(args.a), read(args.b)
    both = sorted(set(A) & set(B))
    print('%s: %d drugs | %s: %d drugs | shared: %d'
          % (args.a_name, len(A), args.b_name, len(B), len(both)))
    if not both:
        return
    ga = sorted({A[d] for d in both})
    gb = sorted({B[d] for d in both})
    ct = collections.Counter((B[d], A[d]) for d in both)
    print()
    w = max(12, max(len(g) for g in gb) + 1)
    print('%-*s' % (w, '%s \\ %s' % (args.b_name, args.a_name)) + ''.join('%11s' % g for g in ga))
    for g2 in gb:
        print('%-*s' % (w, g2) + ''.join('%11d' % ct.get((g2, g1), 0) for g1 in ga))

    if args.a_case and args.b_case:
        n11 = ct.get((args.b_case, args.a_case), 0)
        n10 = ct.get((args.b_case, args.control), 0)
        n01 = ct.get((args.control, args.a_case), 0)
        n00 = ct.get((args.control, args.control), 0)
        print()
        print('2x2 on the case/control cells only: [[%d, %d], [%d, %d]]' % (n11, n10, n01, n00))
        if min(n11, n10, n01, n00) and sps is not None:
            orr = (n11 * n00) / (n10 * n01)
            _, p = sps.fisher_exact([[n11, n10], [n01, n00]])
            print('a %s drug is %.2fx more likely to be %s than a %s drug is  (Fisher p=%.3g)'
                  % (args.b_case, orr, args.a_case, args.control, p))
        else:
            print('a zero cell - no odds ratio')


if __name__ == '__main__':
    main()
