#!/usr/bin/env python

import sys
import argparse
from itertools import groupby
import re
import os

try:
    from Bio import SeqIO
except ImportError:
    msg = """
Could not import the BioPython module which means it is probably
not installed, or at least not available in the PYTHONPATH for this 
particular binary.

If you have conda (recommended) try running:

 $ conda install -c anaconda biopython

or, alternatively with python/pip:

 $ python -m pip install biopython
"""
    sys.stderr.write(msg)
    sys.exit(1)

__author__ = "Joe R. J. Healey"
__version__ = "1.0"
__title__ = "ParseCDHIT"
__license__ = "GPLv3"
__author_email__ = "jrj.healey@gmail.com"


def get_args():
    """Parse command line arguments"""
    desc = """Retrieve cluster sequences from the result of CD-HIT using the .clstr file"""
    epi = """Return the sequences for each cluster defined by a run of CD-HIT.
             A typical CD-HIT run would look like:
             
              $ cd-hit -i seqs.fasta -d 0 -c 0.80 -o clust.fasta
             
             Since this script is reliant on SeqIDs, your fasta SeqIDs must be well formed.
             Ideally this means they are short and unique. The -d 0 flag is required to get
             the SeqIDs as long as possible, thereby minimising the chance that sequences
             get mixed up.
          """

    try:
        parser = argparse.ArgumentParser(description=desc, epilog=epi)
        parser.add_argument('clusterfile', action='store',
                            help='CD-HIT output file (ends in ".clstr").')
        parser.add_argument('fastafile', action='store',
                            help='CD-HIT input file of sequences (a multifasta).')
        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            sys.exit(1)

    except NameError:
        sys.stderr.write("An exception occurred with argument parsing. Check your provided options.")
        sys.exit(1)

    return parser.parse_args()


def main():
    """Call functions and parse the CDHIT file"""

    args = get_args()

    # Add a break/check for duplicate IDs?
    # if idx.contains_duplicates() == True:
    #    sys.exit(1)

    # Get all of the sequence IDs from the cluster file
    with open(args.clusterfile, 'r') as cfh:
        groups = [list(group) for key, group in groupby(cfh, lambda line: line.startswith(">Cluster")) if not key]

    # Load the fastafile ready to be queried
    idx = SeqIO.index(args.fastafile, 'fasta')

    ids = [re.findall(r'>(.*)\.{3}', ''.join(g)) for g in groups]

    for i, cluster in enumerate(ids):
        print("Parsing Cluster %s, with IDs:" % i)
        print(cluster)
        with open('Cluster_'+str(i)+'.fasta', 'w') as ofh:
            for identifier in cluster:
                ofh.write(idx[identifier].format('fasta'))

if __name__ == '__main__':
    main()
