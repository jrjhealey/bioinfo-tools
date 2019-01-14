#!/usr/bin/env python

# Interrogate a sequence by finding out what features are present within a specified range

from __future__ import print_function
import sys
import argparse
from Bio import SeqIO

__author__ = "Joe R. J. Healey"
__version__ = "1.0"
__title__ = ""
__license__ = "GPLv3"
__doc__ = "A parser for the output of MultiGeneBlast"
__author_email__ = "jrj.healey@gmail.com"


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

def get_args():
    """Parse command line arguments"""
    desc = """Probe sequence spans for features."""
    epi = """Interrogate an input sequence for any features that fall within user defined ranges."""

    try:
        parser = argparse.ArgumentParser(description=desc, epilog=epi, prog='get_spans.py')
        parser.add_argument('-v', '--verbose', action='store_true', help='Verbose behaviour (extra information).')
        parser.add_argument('-i', '--infile', action='store', help='Input sequence file.')
        parser.add_argument('-f', '--format', action='store', default='genbank', help='The format of the input file.')
        parser.add_argument('-r', '--range', action='store', help='The range of interest (specified as start:stop).')
        parser.add_argument('-t', '--type', action='store', default='CDS',
                            help='Restrict feature detection to just this type of feature.')

        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            sys.exit(1)

    except NameError:
        sys.stderr.write("An exception occurred with argument parsing. Check your provided options.")
        sys.exit(1)

    return parser.parse_args()


def main():
    args = get_args()
    try:
        rec = SeqIO.read(args.infile, args.format)
    except ValueError:
        sys.stderr.write('Currently this script is only capable of handling a single contiguous file '
                         'not multiple records.')

    start, end = args.range.split(':')
    if args.verbose: print("Start: " + str(start) + ", End: " + str(end), file=sys.stderr)

    desired = set(xrange(int(start),int(end),1))
    for feat in rec.features:
        if feat.type == args.type:
            span = set(xrange(feat.location._start.position, feat.location._end.position))
            if span & desired:
                print(feat)


if __name__ == "__main__":
    main()
