#!/usr/bin/python

# Script to find only exact matches in BLAST tab output.
# Custom columns are supported only if they follow qlen and
# qcovs as columns 13 and 14 (i.e. custom cols from 15 onward).

import sys
import argparse
import traceback
import warnings
import csv
import os
import subprocess
import multiprocessing

# Helper functions
def is_exe(fpath):
	'''Test if PATH entry is an executable file'''
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

def which(program):
	'''Shell `which` emulator to check for BLAST installation in the $PATH'''

	fpath, fname = os.path.split(program)
	if fpath:
		if is_exe(program):
			return program
	else:
		for path in os.environ["PATH"].split(os.pathsep):
			path = path.strip('"')
			exe_file = os.path.join(path, program)
			if is_exe(exe_file):
				return exe_file
	return None


def customBLAST(blastcmd, verbose, basename):
	'''Handle a custom blast command provided directly by the user'''

        if which('blastn') == None or which('blastx') == None or which('blastp') == None or which('tblastx') == None or which('tblastn') == None:
                print("BLAST binaries (blastn/p/x or tblastn/x) don't appear to be present within the $PATH. Check your installation. Exiting.")
                sys.exit(1)

	if 'qlen' not in blastcmd or 'qcovs' not in blastcmd:
       		print('"qlen" as column 13 and "qcov" as column 14 is required for subsequent parsing of exact matches.')
                print('Please reformulate your BLAST command as follows:')
                print('The outfmt should contain: -outfmt "6 std qlen qcovs" (custom columns may follow after).')
                sys.exit(1)
        if '-out' in blastcmd:
		print("The script will generate a temporary blastfile with a predictable filename. DO NOT specify your own in the BLAST command.")
		sys.exit(1)

	resultHandle = "{}.blastout.tmp".format(basename)
	blastcmd = blastcmd + '-out ' + resultHandle

	if verbose == 2:
		print('Running BLAST as follows:' + '\n' + blastcmd)

        BLAST_process = subprocess.Popen(blastcmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        BLAST_process.wait()

	return resultHandle			

# Function that synthesises and runs the BLAST command:
def runBLAST(blastcmd, infile, basename, database, threads, verbose):
	'''Synthesise BLAST commands and make system calls'''

	if which('blastn') == None or which('blastx') == None or which('blastp') == None or which('tblastx') == None or which('tblastn') == None:
		print("BLAST binaries (blastn/p/x or tblastn/x) don't appear to be present within the $PATH. Check your installation. Exiting.")
		sys.exit(1)
	
	resultHandle = "{}.blastout.tmp".format(basename)
	
	blastcmd = 'blastn -query {0} -db {1} -out {2} -perc_identity 100 -num_threads {3} -outfmt "6 std qlen qcovs"'.format(infile, database, resultHandle, threads)
        if verbose == 2:
		print('Running BLAST as follows:' + '\n' + blastcmd)
	
	BLAST_process = subprocess.Popen(blastcmd,shell=True,stdin=subprocess.PIPE,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
	BLAST_process.wait()

	return resultHandle

def parseTabular(row):
	'''Check row of tabular file for match to criteria'''
# 0-based counting of columns
        query = int(0)
        subject = int(1)
        pident = int(2)
        aln_len = int(3)
        mismatch = int(4)
        gaps = int(5)
        qstart = int(6)
        qend = int(7)
        sstart = int(8)
        send = int(9)
        eval = int(10)
        bitscore = int(11)
        qlen = int(12)
        qcovs = int(13)
        # custom columns can be used after qcovs	

	if ( row[gaps]      ==  '0'
	 and row[mismatch]  ==  '0'
         and row[qstart]    ==  '1'
         and row[qlen]      ==  row[qend]
         and row[pident]    ==  '100.000'
         and row[qcovs]     ==  '100'):
		
		return row

def main():
	'''Parse commandline arguments and make function calls for main script functionality.'''
	# Default threads:
	threads = multiprocessing.cpu_count() / 2	

	try:
   		parser = argparse.ArgumentParser(description='This script runs BLAST and returns only exact matches in BLAST tabular outputs. BLAST must be in PATH.')
		parser.add_argument(
			'-i',
			'--infile',
			action='store',
			default=None,
			help='The file to BLAST.')
		parser.add_argument(
			'-o',
			'--outfile',
			action='store',
			help='The filename of the final, filtered, BLAST tabular file. The basename will be used for intermediate files.')
		parser.add_argument(
			'-b',
			'--blastcmd',
			action='store',
			default=None,
			help='A string which corresponds to the custom BLAST command you wish to run (all normal blast options valid). If not specified, a default blastn command will be run.')
		parser.add_argument(
			'-d',
			'--database',
			action='store',
			default=None,
			help='The database file to BLAST sequences against.')
		parser.add_argument(
			'-t',
			'--threads',
			action='store',
			default=threads,
			help='The number of threads BLAST is able to use. Default is half the total machine threads.')
		parser.add_argument(
			'-c',
			'--compliant',
			action='store_true',
			help='If provided, the script will additionally output a standard "outfmt 6" file of exact matches for compliance with later steps.')
		parser.add_argument(
			'-p',
			'--parseonly',
			action='store',
			default=None,
			help='If an existing BLAST tab file with correct columns is passed to --parseonly, BLAST will not be rerun, instead exact matches from the provided file will be found.')
		parser.add_argument(
			'-v',
			'--verbose',
			type=int,
			choices=[0,1,2],
			default=0,
			help='Degree of verbosity: Specify 0 [default] for no output. 1 to print the filtered result. 2 to print filtered result and intermediate messages.')
		args = parser.parse_args()
	except:
                print("An exception occured with argument parsing. Check your provided options.")
		traceback.print_exc()

	infile = args.infile	
	if args.outfile is None:
		args.outfile = os.path.basename(os.path.splitext(infile)[0])
	split     = os.path.splitext(args.outfile)
	basename  = os.path.basename(split[0])
	blastcmd  = args.blastcmd
	database  = args.database
	parseonly = args.parseonly
	verbose   = args.verbose

# Main code:

	# If not just parsing an existing file:
	if parseonly is None:
		resultHandle = runBLAST(blastcmd, infile, basename, database, threads, verbose)
		with open(resultHandle, 'r') as tabfile:
			tsvin = csv.reader(tabfile, delimiter='\t')
			exacts = []
			for row in tsvin:
				match = parseTabular(row)
				exacts.append(match)

	elif parseonly is not None:
		with open(parseonly, 'r') as tabfile:
                	tsvin = csv.reader(tabfile, delimiter='\t')
                	exacts = []
			for row in tsvin:
				match = parseTabular(row)
				exacts.append(match)

# 

	with open(basename + '_exacts.tsv', 'w') as ofh:
		for each in exacts:
			if each is not None:
				ofh.write('\t'.join(each))
#	if compliant == True:
#		with open(basename + '_exacts_compliant.tsv', 'w') as cfh:
#			for each in exacts:
#				if each is not None:
#					cfh.write('\t'.join(each[0:12])

	if verbose != 0:
		print("\n")
		for each in exacts:
			if each is not None:
				print('\t'.join(each))

if __name__ == '__main__':
	main()
