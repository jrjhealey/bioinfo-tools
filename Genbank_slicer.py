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
import os
		
def convert(basename, genbank):
	'''Convert the provided genbank to a fasta to BLAST.'''
	
	refFasta = "{}.fasta.tmp".format(basename)
	SeqIO.convert(genbank, 'genbank', refFasta, 'fasta')

	return refFasta

def runBlast(basename, refFasta, fasta):
	'''Synthesise BLAST commands and make system calls'''
	
	resultHandle = "{}.blastout.tmp".format(basename)
	blastdb_cmd = 'makeblastdb -in {0} -dbtype nucl -title temp_blastdb'.format(refFasta)
	blastn_cmd = 'blastn -query {0} -strand both -task blastn -db {1} -perc_identity 100 -outfmt 6 -out {2} -max_target_seqs 1'.format(fasta, refFasta, resultHandle)
	print("Constructing BLAST Database: " + '\n' + blastdb_cmd)
	print("BLASTing: " + '\n' + blastn_cmd)
	DB_process = subprocess.Popen(blastdb_cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	DB_process.wait()
	BLAST_process = subprocess.Popen(blastn_cmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	BLAST_process.wait()

	return resultHandle
	     
def getIndices(resultHandle):
	'''If not provided directly by the user, this function retrieves the best BLAST hit's indices.'''
	
	blast_result = SearchIO.read(resultHandle, 'blast-tab')
	print(blast_result[0][0])
	start = blast_result[0][0].hit_start
	end = blast_result[0][0].hit_end
		
	return start, end

def slice(start, end, genbank):
	'''Subset the provided genbank to return the sub record.'''
	
	seqObj = SeqIO.read(genbank, 'genbank')
	subRecord = seqObj[start:end]

	return subRecord


def main():
###################################################################################################
# Parse command line arguments
	try:
   		parser = argparse.ArgumentParser(
   			description='This script slices entries such as genes or operons out of a genbank, subsetting them as their own file.')
		parser.add_argument(
			'-g',
			'--genbank',
			action='store',
			required=True,
			help='The genbank file you wish to subset.')
		parser.add_argument(
			'-o',
			'--outfile',
			action='store',
			help='If specifed, the script will write a file, otherwise redirect STDOUT for pipelines.')
		parser.add_argument(
			'-s',
			'--start',
			type=int,
			help='Integer base index to slice from.')
		parser.add_argument(
			'-e',
			'--end',
			type=int,
			help='Integer base index to slice to.')
		parser.add_argument(
			'-b',
			'--blastmode',
			action='store_true',
			default=False,
			help='(False if not specified- If flag is provided, you can provide an input fasta, and BLAST will be called to find your sequence indices in the genbank')
		parser.add_argument(
			'-f',
			'--fasta',
			action='store',
			help='The operon fasta to pull annotations from the provided genbank.')
		parser.add_argument(
			'-t',
			'--no_tidy',
			action='store_false',
			default=True,
			help='Tell the script whether or not to remove the temporary files it generated during processing. On by default. WARNING: removes files based on the "tmp" string.')
		
		args = parser.parse_args()

	except:
		print "An exception occured with argument parsing. Check your provided options."
		traceback.print_exc()
	
	genbank   =  args.genbank
	fasta     =  args.fasta
	split     =  os.path.splitext(args.genbank)
	basename  =  os.path.basename(split[0])
	start     =  args.start
	end       =  args.end
	blastMode =  args.blastmode
	outfile   =  args.outfile
	fasta     =  args.fasta
	tidy      =  args.no_tidy


# Main code:
	if blastMode is not False and fasta is not None:
		refFasta = convert(basename,genbank)
		resultHandle = runBlast(basename, refFasta, fasta)
		start, end = getIndices(resultHandle)
	else:
		if fasta is None:
			print("No fasta was provided so BLAST mode cannot be used.")
      
	if start is None or end is None:
        	print('No slice indices have been specified or retrieved from blastout')
                sys.exit(1)


	subRecord = slice(start, end, genbank)
	
	
	if outfile is not None:
		SeqIO.write(subRecord, outfile, "genbank")
	else:
		print(subRecord.format('genbank'))

	if tidy is True:
		subprocess.Popen("rm -v ./*tmp*",shell=True)

if __name__ == "__main__":
	main()

