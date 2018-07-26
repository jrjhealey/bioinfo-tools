#/usr/bin/env python

"""
Get the distances between nodes in a phylogenetic tree.

Requires dendropy.
"""
from __future__ import print_function
import argparse, sys
try:
   import dendropy
except ImportError:
    msg = """
The dendropy import failed, the module doesn't appear to be installed
(at least in the PYTHONPATH for this python binary").
Try running:
 $ python -m pip install -U dendropy
or 
 $ conda install -c etetoolkit ete3 ete_toolchain
"""
    print(msg)
    sys.exit(1)

def get_args():
    """Parse command line arguments."""
    try:
        parser = argparse.ArgumentParser(description="Get the distances between nodes in a tree.",
                                     formatter_class=argparse.RawTextHelpFormatter)
        parser.add_argument('-i','--input', action='store', required=True, help='Input tree file.')
        parser.add_argument('-s', '--schema', action='store', required=True, help='Input tree format.')
        parser.add_argument('-m', '--mode', action='store', choices=['all', 'max'], default='max',
                            help='What distances to return: Largest = \'max\', or All = \'all\'.')
        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            sys.exit(1)
    except:
        print("An exception occurred with argument parsing. Check your provided options.", file=sys.stderr)
        sys.exit(1)

    return parser.parse_args()


def err(string):
    """Print errors to stderr"""
    print("%s" % string, file=sys.stderr)


def main():
    """Calculate the largest and smallest distances in a phylogenetic tree via Dendropy."""
    args = get_args()

    # Run dendropy functions
    tree = dendropy.Tree.get(path=args.input, schema=args.schema)
    pdm = tree.phylogenetic_distance_matrix()
    if args.mode is 'max':
        max1, max2 = pdm.max_pairwise_distance_taxa()
        print("Distance between '%s' and '%s': %s" % (max1.label, max2.label, pdm(max1, max2)))
    else:
        for i, t1 in enumerate(tree.taxon_namespace[:-1]):
            for t2 in tree.taxon_namespace[i+1:]:
                print("Distance between '%s' and '%s': %s" % (t1.label, t2.label, pdm(t1, t2)))


if __name__ == "__main__":
    main()
