import omnipath as op
from omnipath_client.utils import map_name
import argparse
import os

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out_dir", default="./")

    args = ap.parse_args()


    interactions = op.interactions.OmniPath.get(genesymbols=True)
    enzsub = op.requests.Enzsub.get(genesymbols=True)
    collectri = op.interactions.CollecTRI.get(
        organism="human",
        genesymbols=True,
        loops=True
    )

    interactions.to_csv(os.path.join(args.out_dir, "interactions.tsv"), sep = '\t')
    enzsub.to_csv(os.path.join(args.out_dir, "enzsub.tsv"), sep = '\t')
    collectri.to_csv(os.path.join(args.out_dir, "collectri.tsv"), sep = '\t')
if __name__ == "__main__":
    main()
