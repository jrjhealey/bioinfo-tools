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

#TODO:
# Add some regex support to the Hamming function to permit any char?
# Add control/logic for case-sensitivity

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
			    type=str,
			    action='store',
			    help='Pass the first string to be compared directly as text.')
	parser.add_argument('-B',
			    '--stringB',
			    type=str,
			    action='store',
			    help='Pass the second string to be compared directly as text.')
	parser.add_argument('--average',
			    action='store_false',
			    help='Average the distance between all sequences (for MSAs).')
        parser.add_argument('-v',
                            '--verbose',
                            action='store_true',
                            help='Prints additional progress messages.')

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
    """Compute Hamming distance from a provided MSA or pair of strings."""

    args = parseArgs()

    if args.alignment is not None:
        if args.verbose: print("Alignment found, returning all pairwise distances")
        from Bio import AlignIO
        msa = AlignIO.read(args.alignment, args.format)
        dists = []
        seq1_list = []
        seq2_list = []
        for i in xrange(len(msa)):
            for j in xrange(i+1, len(msa)):
                dists.append(hamming_distance(str(msa[i].seq), str(msa[j].seq)))
                seq1_list.append(str(msa[i].seq))
                seq2_list.append(str(msa[j].seq))
        print("Hamming distances:")
        for dist, seq1, seq2 in zip(dists, seq1_list, seq2_list):
            print(str(dist) + '\t' + seq1 + '\t' + seq2)

    else:
        if args.verbose: print("No MSA found, comparing strings instead.")
	if args.stringA and args.stringB:
	    ham_pair = hamming_distance(args.stringA, args.stringB)
	    print(str(ham_pair) + '\t' + args.stringA + '\t' + args.stringB)
        else:
	    print("Strings A and B are undefined, ensure both are provided with -A/--stringA and -B/--stringB")


if __name__ == '__main__':
    main()
