#!/usr/bin/python

# This script is designed to take a genbank file and 'slice out'/'subset'
# regions (genes/operons etc.) and produce a separate file. This can be
# done explicitly by telling the script which base sites to use, or can
# 'decide' for itself by blasting a fasta of the sequence you're inter-
# ed in against the Genbank you want to slice a record out of.

# Note, the script (obviously) does not preseve the index number of the
# bases from the original 

# Based upon the tutorial at:
# http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc44

# This script depends on BLAST and having the makeblastdb functionality
# available if BLAST_MODE is active. It also depends on Biopython.

# Set up and handle arguments:

from Bio import SeqIO
from Bio import SearchIO
import subprocess
import sys
import argparse
import traceback
import warnings


def main():
###################################################################################################
# Parse command line arguments
	try:
   		parser = argparse.ArgumentParser(
   			description='This script slices entries such as genes or operons out of a genbank, subsetting them as their own file.')
		parser.add_argument(
			'-g',\
			'--genbank',\
			action='store',\
			required=True,\
			help='The genbank file you wish to subset.')
		parser.add_argument(
			'-o',\
			'--outfile',\
			action='store',\
			help='If specifed, the script will write a file, otherwise redirect STDOUT for pipelines.')
		parser.add_argument(
			'-s',\
			'--start',\
			type=int,\
			help='Integer base index to slice from.')
		parser.add_argument(
			'-e',\
			'--end',\
			type=int,\
			help='Integer base index to slice to.')
		parser.add_argument(
			'-b',
			'--blast_mode',\
			action='store_true',\
			default=False,\
			help='(True|False - If true, you can provide an input fasta, and BLAST will be called to find your sequence indices in the genbank')
		parser.add_argument(
			'-f',\
			'--fasta',\
			action='store',\
			help='The operon fasta to pull annotations from the provided genbank.')
		args = parser.parse_args()
	except:
		print "An exception occured with argument parsing. Check your provided options."
		traceback.print_exc()



# Main code:
	if args.blast_mode is not False:
		split = os.path.splitext(args.genbank)
		basename = os.path.basename(split[0])
		ref_fasta = "{}.fasta.tmp".format(basename)
		result_handle = "{}.blastout.tmp".format(basename)
                SeqIO.convert(args.genbank, 'genbank', ref_fasta, 'fasta')
		blastdb_cmd = 'makeblastdb -in {0} -dbtype nucl -title temp_blastdb'.format(ref_fasta)
		blastn_cmd = 'blastn -query {0} -strand both -task blastn -db {1} -perc_identity 100 -outfmt 6 -out {2} -max_target_seqs 1'.format(args.fasta, ref_fasta, result_handle)
		print("Constructing BLAST Database: " + blastdb_cmd)
		print("BLASTing: " + blastn_cmd)
		subprocess.Popen(blastdb_cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		subprocess.Popen(blastn_cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
		blast_result = SearchIO.read(result_handle, 'blast-tab')
		print(blast_result[0][0])
		args.start = blast_result[0][0].hit_start
		args.end = blast_result[0][0].hit_end
		seq_obj = SeqIO.read(args.genbank, 'genbank')
		sub_record = seq_obj[args.start:args.end]
		if args.outfile is not None:
            	  	SeqIO.write(sub_record, args.outfile, "genbank")
		else:
			print(sub_record.format('genbank'))
		
        else:
                if args.start is None or args.end is None:
                        print('No slice indices have been specified or retrieved from blastout')
                        sys.exit(0)
                else:
                        seq_obj = SeqIO.read(args.genbank, 'genbank')
			sub_record = seq_obj[args.start:args.end]
			if args.outfile is not None:
                        	SeqIO.write(sub_record, args.outfile, "genbank")
			else:
				print(sub_record.format('genbank'))


if __name__ == "__main__":
	main()

