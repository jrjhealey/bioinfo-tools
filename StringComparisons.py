#/usr/bin/env python

# Calculate similarity scores between sets of sequences via different metrics.
# Requires BioPython as the only non-standard module.

import os
import sys
import argparse
import re
import math
from functools import partial
from collections import Counter
from Bio import AlignIO

__author__ = "Joe R. J. Healey"
__version__ = "1.1"
__title__ = "StringComparisons.py"
__license__ = "GPLv3"
__author_email__ = "jrj.healey@gmail.com"

#TODO:
# Add regex/fuzzy matching to the distance functions?
# Add control/logic for case-sensitivity

def get_args():
    """Parse command line arguments"""
    desc = 'Perform various string comparisons using different metrics.'
    epi = ('This is a little script to perform various string comparisons '
          'between elements of a sequence alignment')

    try:
        parser = argparse.ArgumentParser(description=desc, epilog=epi, prog='StringComparisons.py')

        parser.add_argument('-a',
                            '--alignment',
                            action='store',
                            help='The multiple sequence alignment (MSA) in any of the formats supported by Biopython\'s AlignIO.')
	parser.add_argument('-f',
			    '--format',
			    action='store',
			    default='fasta',
			    help='The format of the sequence alignment, if not the default = FASTA.')
        parser.add_argument('-m',
                            '--method',
                            action='store',
                            choices=['hamming','cosine','levenshtein'],
                            default='levenshtein',
                            metavar="METHOD",
                            help='What type of string comparison measure to return {hamming|cosine|levenshtein} [default = levenshtein]'
                                 'If cosine is chosen, an optional kmer length is needed via -k|--kmer')
        parser.add_argument('-k',
                            '--kmer',
                            type=int,
                            default=5,
                            help='Kmer length for use with method = "cosine".')
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
        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            sys.exit(1)

    except NameError:
        sys.stderr.write("An exception occurred with argument parsing. Check your provided options.")
        sys.exit(1)

    return parser.parse_args()


def hamming_distance(s1, s2):
     """Return the Hamming distance between equal-length sequences"""
     if len(s1) != len(s2):
         raise ValueError("Undefined for sequences of unequal length")
     return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def levenshtein_distance(s1, s2):
    if len(s1) < len(s2):
        return levenshtein_distance(s2, s1)

    # len(s1) >= len(s2)
    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    return previous_row[-1]

def cosine_distance(s1, s2, k):
    """Compute the cosine difference of the strings as kmer vectors"""
    import re, math
    from collections import Counter

    # Convert DNA "sentence" to 'word' (kmer) vector
    kmers = re.compile("(?=(\w{%s}))" % k)
    vec1 = Counter(kmers.findall(s1))
    vec2 = Counter(kmers.findall(s2))

    intersection = set(vec1.keys()) & set(vec2.keys())
    numerator = sum([vec1[x] * vec2[x] for x in intersection])

    sum1 = sum([vec1[x]**2 for x in vec1.keys()])
    sum2 = sum([vec2[x]**2 for x in vec2.keys()])
    denominator = math.sqrt(sum1) * math.sqrt(sum2)
    if not denominator:
        return 0.0
    else:
        return float(numerator) / denominator

def apply_method(method, s1, s2, k):
    """Case switch for the selected method"""
    return {
            'hamming': partial(hamming_distance, s1, s2),
            'cosine': partial(cosine_distance, s1, s2, k),
            'levenshtein': partial(levenshtein_distance, s1, s2)
           }[method]()

def main():
    """Compute distances from a provided MSA or pair of strings."""

    args = get_args()

    if args.alignment is not None:
        if args.verbose: print('Alignment found, returning all pairwise distances.')
        msa = AlignIO.read(args.alignment, args.format)
        dists = []
        seq1_list = []
        seq2_list = []
        for i in range(len(msa)):
            for j in range(i+1, len(msa)):
                dists.append(apply_method(args.method, str(msa[i].seq), str(msa[j].seq), k=args.kmer))
                seq1_list.append(str(msa[i].seq))
                seq2_list.append(str(msa[j].seq))
        for dist, seq1, seq2 in zip(dists, seq1_list, seq2_list):
            print(str(dist) + '\t' + seq1 + '\t' + seq2)
        if args.average:
            print("Average pairwise similarity between MSA: {}".format(sum(dists)/len(dists)))
    else:
        if args.verbose: print("No MSA found, comparing strings instead.")
	if args.stringA and args.stringB:
	    result = apply_method(args.method, args.stringA, args.stringB, args.kmer)
	    print(str(result) + '\t' + args.stringA + '\t' + args.stringB)
        else:
	    sys.stderr.write("Strings A and B are undefined, ensure both are provided with -A/--stringA and -B/--stringB")


if __name__ == '__main__':
    main()
