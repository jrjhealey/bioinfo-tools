#!/home/wms_joe/bin/miniconda2/bin/python
"""
This script will run and parse in HHpred results.
The results will be added to an easily readable tsv.
"""

import os
import subprocess
import sys
import argparse
import traceback
import warnings
import pandas as pd
from io import StringIO

# HHSuite output parser/template

template = \
u"""
---|------|------------------------|----|-------|-------|------|-----|----|---------|--------------|
 No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
"""

def parseTable(hhresult_file):
    	pattern = StringIO(template).readlines()[1]
    	colBreaks = [i for i, ch in enumerate(pattern) if ch == '|']
    	widths = [j-i for i, j in zip( ([0]+colBreaks)[:-1], colBreaks ) ]

    	hhtable = pd.read_fwf(hhresult_file, skiprows=8, nrows=10, header=0, widths = widths)

    	top_hit = hhtable.loc[0,'Hit']
    	top_prob = hhtable.loc[0,'Prob']
    	top_eval = hhtable.loc[0,'E-value']
    	top_pval = hhtable.loc[0,'P-value']
    	top_score = hhtable.loc[0,'Score']

    	return top_hit, top_prob, top_eval, top_pval, top_score

def parseHeader(hhresult_file):
	with open(hhresult_file, 'r') as rfh:
		for line in rfh:
			if line.startswith('>'):
				header = line.rstrip('\n')

	return header
	
##### Main code begins #####



# Step 1, parse the output of HHpred to get the nearest homolog.

def main():
	def_cpus=1
	def_dir=os.getcwd()
	try:
   		parser = argparse.ArgumentParser(description='This script will run HHpred and parse the output in to a TSV file.')

		parser.add_argument(
			'-d',
			'--database',
			action='store',
			default=None,
			help='You can specify a different HMM database (filepath) to use if you offer this parameter. Otherwise it defaults to PDB.')

		parser.add_argument(
			'-f',
			'--fasta',
			action='store',
			help='The fasta amino acid sequence that corresponds to the simulated structure and the sequence you wish to query.')

		parser.add_argument(
			'-t',
			'--threads',
			action='store',
			default=def_cpus,
			help='The number of threads that HHsearch can execute on.')

		parser.add_argument(
			'-o',
			'--outdir',
			action='store',
			default=def_dir,
			help="Directory for files to be written to. Default is the current working directory.")
		args = parser.parse_args()

	except:
		print("An exception occured with argument parsing. Check your provided options.")
		traceback.print_exc()

  # Acquire HHPred Database
	if args.database is None:
		db_cmd = 'find ~/Applications/HHSuite/databases/pdb70 -type f -name "pdb70_hhm.ffdata"'
        	db_path = subprocess.Popen(
                	db_cmd,
                	shell=True,
                	stdin=subprocess.PIPE,
                	stdout=subprocess.PIPE,
                	stderr=subprocess.PIPE)

		stdout, stderr = db_path.communicate()
        	filelist = stdout.decode().split()
		args.database = stdout
	else:
		print("No default database found and you haven't provided one directly.")

	if args.fasta is not None:
		fastaDir = os.path.splitext(args.fasta)
		basename = os.path.basename(fastaDir[0])
	else:
		print('No input fasta was provided. Check your input parameters')
		sys.exit(1)

	# Run the HHsearch
	if args.fasta is not None and args.database is not None:
		print("\n")
		print("Running " + args.fasta + " in HHsearch, against " + args.database + " on " + str(args.threads) + " threads.")
		hhresult_file = '{0}.hhr'.format(os.path.join(os.path.abspath(args.outdir),basename))
		search_cmd = 'hhsearch -dbstrlen 50 -B 1 -b 1 -p 60 -Z 1 -E 1E-03 -nocons -nopred -nodssp -cpu {0} -i {1} -d {2} -o {3}'.format(args.threads,args.fasta,args.database,hhresult_file)
		print("\n")
		print("Executing HHsearch with the command:")
		print(search_cmd)
		hh_process = subprocess.Popen(
		 			search_cmd,
					shell=True,
					stdin=subprocess.PIPE,
					stderr=subprocess.PIPE)
		hh_process.wait()

	else:
		print("No fasta or database has been detected. Check your input parameters.")
		sys.exit(1)

	# Parse HHpred output and acquire the best hit

        top_hit, top_prob, top_eval, top_pval, top_score = parseTable(hhresult_file)
	header = parseHeader(hhresult_file)	

        print("Your best hit: (PDB ID | Probability | E-Value | P-Value | Score)")
        print("\t" +  str(top_hit) + "\t" + str(top_prob) + "\t" + str(top_eval) + "\t" + str(top_pval) + "\t" + str(top_score) )
	print(header)

	print("Full results are found in: " + hhresult_file)

	outfile = '{}.tsv'.format(os.path.join(args.outdir, basename))
	
	with open(outfile, 'w') as ofh:
		ofh.write(basename + "\t" \
                          + str(top_hit) + "\t" \
                          + str(header) + "\t" \
                          + str(top_prob) + "\t" \
                          + str(top_eval) + "\t" \
                          + str(top_pval) + "\t" \
                          + str(top_score) + "\n")

	print('All done.')
	sys.exit(0)

if __name__ == '__main__' :
	main()
