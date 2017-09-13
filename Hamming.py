# script to calculate the Hamming distance between elements of a MSA


import os
import sys
import warnings
import traceback

__author__ = "Joe R. J. Healey"
__version__ = "1.0.0"
__title__ = "Hamming.py"
__license__ = "GPLv3"
__author_email__ = "J.R.J.Healey@warwick.ac.uk"


def parseArgs():
    """Parse command line arguments"""

    import argparse

    try:
        parser = argparse.ArgumentParser(
            description='Calculate the Hamming Distance between sets of alignments.')

        parser.add_argument('-a',
                            '--alignment',
                            action='store',
                            help='The multiple sequence alignment (MSA) in any of the formats supported by Biopython\'s AlignIO.')
	parser.add_argument('-f',
			    '--format',
			    action='store',
			    default='fasta',
			    help='The format of the sequence alignment, if not the default = FASTA.')
	parser.add_argument('-A',
			    '--stringA',
			    action='store',
			    help='Pass the first string to be compared directly as text.')
	parser.add_argument('-B',
			    '--stringB',
			    action='store',
			    help='Pass the second string to be compared directly as text.')
	parser.add_argument('--average',
			    action='store_false',
			    help='Average the distance between all sequences (for MSAs).')

    except:
        print "An exception occurred with argument parsing. Check your provided options."
        traceback.print_exc()

    return parser.parse_args()


def hamming_distance(s1, s2):
     """Return the Hamming distance between equal-length sequences"""
     if len(s1) != len(s2):
         raise ValueError("Undefined for sequences of unequal length")
     return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def main():
    """Compute Shannon Entropy from a provided MSA."""

    args = parseArgs()

    if args.alignment is not None:
    	from Bio import AlignIO
    	msa = AlignIO.read(args.alignment, args.format)
	for a, i in enumerate(msa):
		for b, j in enumerate(msa):
			

	else:
		# Report pairwise

# If not alignment: do pairwise, return un-averaged



if __name__ = '__main__':
    main()
