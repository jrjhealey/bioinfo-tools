"""
Calculate dumb and brute force consensus sequences from MSAs.
Currently requires that the MSA has no entire gap columns.
"""

# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

import os
import sys
import traceback
import warnings
from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

__author__ = "Joe R. J. Healey"
__version__ = "1.0.0"
__title__ = "Consensus.py"
__license__ = "GPLv3"
__author_email__ = "J.R.J.Healey@warwick.ac.uk"


def parse_args():
	"""Parse commandline arguments"""
	import argparse
	from argparse import RawTextHelpFormatter
	from textwrap import dedent
	try:
		parser = argparse.ArgumentParser(
				description='Calculate dumb and brute force consensus sequences from Multiple Sequence Alignments.',
				usage='python Consensus.py [options] task alignment\nRun with all defaults: python Consensus.py alignment\n\npython Consensus.py -h|--help for full options',
				formatter_class=RawTextHelpFormatter)
		parser.add_argument('task',
							action='store',
							default='forced',
							choices=['dumb','forced','pssm'],
							const='forced',
							nargs='?',
							help='Which task to perform.\nForced will randomly resolve ambiguous positions from the most likely choices [Default].')
		parser.add_argument('alignment',
							action='store',
							help='The MSA to analyse.')
		parser.add_argument('-f',
							'--format',
							action='store',
							default='None',
							help='Alignment file format.\nScript will attempt to guess from the extension but may fail.')
		parser.add_argument('-v',
							'--verbose',
							action='store_true',
							help='Print additional messages to screen. [False].')
		parser.add_argument('-c',
							'--cout',
							action='store',
							default=None,
							help='Output file for consensus sequence.\nPrints to screen in fasta format by default.')
		# parser.add_argument('-m',
		# 					'--matrix',
		# 					action='store',
		# 					default=None,
		# 					choices=['dumb', 'forced', 'N'],
		# 					help=dedent('''\
		# 					Output the Position Specific Score Matrix of the consensus.
		# 					Choose a representative sequence for the PSSM axis from:
		# 					 "dumb" consensus   -  include ambiguous characters.
		# 					 "forced" consensus -  ambiguous characters are resolved randomly.
		# 					 "N"                -  some integer specifying a sequence in the alignment.
		# 					'''))
		parser.add_argument('-p',
							'--pout',
							action='store',
							default=None,
							help='Output file to store the Position Specific Score Matrix if -m|--matrix was given. Else prints to screen.')
		# Arguments to BioPython's dumb_consensus
		parser.add_argument('-t',
							'--threshold',
							action='store',
							type=int,
							default=0.7,
							help='Frequency threshold for inclusion of residue in to consensus, passes through to the dumb_consensus method.')
		parser.add_argument('-a',
							'--ambiguous',
							action='store',
							default='X',
							help='The ambiguous character used for dumb consensuses.')

	except NameError:
		print "An exception occured with argument parsing. Check your provided options."
		traceback.print_exc()
		sys.exit(1)

	return parser.parse_args()



def guess_ext(args):
	"""If no input ext was specified, guess it from the extension name or emit a warning"""
	extension = os.path.splitext(args.alignment)[1]
	if args.verbose is True: print("Extension is " + extension)
	# Figure out what extension to return
	if extension in (".clust", ".clustal",".aln"):
		ext = "clustal"
	elif extension in (".emb", ".emboss"):
		ext = "emboss"
	elif extension in (".fa", ".fasta", ".fas", ".fna", ".faa", ".afasta"):
		ext = "fasta"
	elif extension in (".phy", ".phylip"):
		ext = "phylip"
	elif extension in (".nexus", ".paup"):
		ext = "nexus"
	else:
		print("Couldn't determine the file format from the extension. Reattempt with the -f|--format option specified.")
		print("Acceptable formats are those supported by BioPython AlignIO. Comprehensive list at http://biopython.org/wiki/AlignIO. Exiting.")
		sys.exit(1)

	if args.verbose is True: print("Your file looks like:" + ext)

	return ext

def dumb_cons(args, msa_summary):
	"""Compute a dumb consensus sequence using BioPython"""
	return SeqRecord(msa_summary.dumb_consensus(threshold=args.threshold,
												ambiguous=args.ambiguous),
					 id=os.path.basename(args.alignment) + ' consensus',
					 description='',
					 name='')


def brute_force_cons(args, msa_summary):
	"""Take a consensus sequence similar to Biopython's dumb one,
	but randomly resolve equal likelihood residues. This is a majority-rule
	consensus, not taking in to consideration thresholds, unlike the dumb_consensus."""

	import random


	# Prepare the string to be the consensus
	consensus = ''
	# Iterate the full MSA length column by column
	for i in xrange(msa_summary.alignment.get_alignment_length()):
		possibles = enumerate_string(msa_summary.get_column(i))
		if len(possibles) == 1:
			consensus += possibles[0]
		else:
			consensus += random.choice(possibles)

	consensus_record = SeqRecord(Seq(consensus),
								 id=os.path.basename(args.alignment) + ' consensus',
								 description='',
								 name='')

	return consensus_record

def enumerate_string(string):
	"""Returns the most common characters of a string. Multiple characters are returned if there are equally frequent characters."""

	from collections import Counter
	# Get counts of each letter in the string
	counts = Counter(string)

	keys = []
	# Iterate each key value pair in the counts dict in case there are tied values
	for key, value in counts.iteritems():
		# collect all keys which occur equally many times
		if value == max(counts.values()):
			keys.append(key)

	return keys



def main():
	"""Call functions to perform MSA analyses."""

	args = parse_args()

	if args.format is "None": args.format = guess_ext(args)
	msa = AlignIO.read(args.alignment, args.format)
	msa_summary = AlignInfo.SummaryInfo(msa)

	# Switch for which action is performed

	if args.task == "pssm":
		pssm = msa_summary.pos_specific_score_matrix()
		if args.pout is not None:
			import csv
			with open(args.pout, 'wb') as f:
				w = csv.DictWriter(f, pssm[0].keys())
				w.writerow(dict((fn, fn) for fn in pssm[0].keys()))
				w.writerows(pssm)
		else:
			print(pssm)
		sys.exit(0)

	elif args.task == "dumb":
		consensus = dumb_cons(args, msa_summary)

	elif args.task == "forced":
		consensus = brute_force_cons(args, msa_summary)

	if args.cout is not None:
		SeqIO.write(consensus, args.cout, 'fasta')
	else:
		print(consensus.format('fasta'))

if __name__ == '__main__':
	main()
