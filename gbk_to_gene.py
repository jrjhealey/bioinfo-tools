# Script to get a feature from a Genbank file

import traceback
import os
import warnings

#TODO
# Set script up to take arguments of the form locus_tag:name or product:string
# so that multiple things can be searched for.

__author__ = "Joe R. J. Healey"
__version__ = "1.0.0"
__title__ = "fetch.py"
__license__ = "GPLv3"
__author_email__ = "J.R.J.Healey@warwick.ac.uk"

def parseArgs():
	"""Parse command line arguments"""

	import argparse

	try:
		parser = argparse.ArgumentParser(description="Retrieve features from a genbank by name.')
		parser.add_argument('feature',
                            	    action='store',
                            	    help='The type and string of the feature you want colon separated, e.g. "locus_tag:EC_1000" or "product:dnaA".')
		parser.add_argument('--from',
                            	    action='store',
                            	    help='The genbank you wish to pull the feature from.')
                parser.add_argument('-t',
                		   '--type',
                		   choices=('nuc','prot'),
                		   action='store',
                		   default='nuc',
                		   help='What type of feature to return, the DNA or translated protein sequence. Nucleotide by default.')
	except:
		print "An exception occurred with argument parsing. Check your provided options."
		traceback.print_exc()

    return parser.parse_args()


def fetchFeature(genbank, feature, type):
	"""Pull a feature from a Genbank, by name"""
	
	from Bio import SeqIO
	seq_record = SeqIO.read(genbank, 'genbank')
	for seq_feature in seq_record.features:
		if seq_feature.type == 'CDS':
			if seq_feature.qualifiers['gene'][0] == feature or seq_feature
			assert len(seq_feature.qualifiers['translation'])==1
                        seq_feature.qualifiers['gene'][0],
                                seq_record.name,
                                seq_feature.qualifiers['translation'][0]))


def main():
	args = parseArgs()
	genbank = args.genbank
	feature = args.feature
	type = args.type

	fetchFeature(genbank, feature)

if __name__ == "__main__":
	main
