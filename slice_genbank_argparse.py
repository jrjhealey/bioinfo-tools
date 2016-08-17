#!/usr/bin/python

# This script is designed to take a genbank file and 'slice out'/'subset'
# regions (genes/operons etc.) and produce a separate file.

# Based upon the tutorial at:
# http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc44

# This script depends on BLAST and having the makeblastdb functionality
# available if BLAST_MODE is active. It also depends on Biopython.

# Set up and handle arguments:

from Bio import SeqIO
from Bio import SearchIO
from Bio.Blast.Applications import NcbiblastnCommandline
from subprocess import call
import os
import sys
import argparse
import traceback

###

def main():
	if blast_mode is not None:	# Test for blast-mode
		conv_fasta = SeqIO.convert(genbank, 'genbank', genbank[:-3] +'tmp', 'fasta')
		blastdb = "makeblastdb -in tmp.fasta -dbtype 'nucl' -title genbank[:-3] +'db'"
		os.system(blastdb)
		os.system(blast_cmd) # blastn -query query -strand both -task blastn -db ./genbank
		# Blast input sequences against the database to get the slice indexes
		# Converting the genbank
		print('BLAST option provided. Executing:')
		print(blastn_cline)
		# Runs the blast
		blastn_cline
		stdout, stderr = blastn_cline()
		# Read in the blastout produced above and parse Name, ID and indices.
		dbresult = SearchIO.read({record}.blastout, 'blast-tab')
		index_start = dbresult.colx
		index_end = dbresult.coly
	else:	
		SeqIO.read(genbank, 'genbank')
	if index_start is None or index_end is  None:
		print('No slice indices have been specified or retrieved from blastout')
		sys.exit(0)
	else:
		print("Slicing " + filename + " from " + str(index_start) + " to " + str(index_end))
		sub_record = genbank[start:end]
		SeqIO.write(sub_record, outfile, "genbank")

if __name__ == '__main__':
	try:
   		parser = argparse.ArgumentParser(
   			description='This script slices entries such as genes or operons out of a genbank, subsetting them as their own file.')
		parser.add_argument(
			'-g', '--genbank', action='store', required=True,\
			help='The genbank file you wish to subset.')
		parser.add_argument(
			'-o', '--outfile', action='store',\
			help='If specifed, the script will write a file, otherwise redirect STDOUT for pipelines.')
		parser.add_argument(
			'-s', '--start', type=int,\
			help='Integer base index to slice from.')
		parser.add_argument(
			'-e', '--end', type=int,\
			help='Integer base index to slice to.')
		parser.add_argument(
			'-b', '--blast_mode', action='store_true', default=False,\
			help='If this switch is toggled, you can provide an input fasta, and BLAST will be called to find your sequence indices in the genbank')
		parser.add_argument(
			'-f', '--fasta', action='store',\
			help='The operon fasta to pull annotations from the provided genbank.')
		args = parser.parse_args()
		genbank = SeqIO.read(args.genbank, "genbank")
		filename = args.genbank
		outfile = args.outfile
		index_start = args.start
		index_end = args.end
		query = args.fasta
		blast_mode = args.blast_mode
		print args
	except:
		print "Issues with provided args."
		traceback.pring_exc()
