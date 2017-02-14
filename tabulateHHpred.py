
# This scripts takes the overly verbose HHpred outputs
# and makes tabulated outputs

import os
import subprocess
import sys
import argparse
import traceback
import warnings
import pandas as pd
from io import StringIO

# Template of HHpred's verbose tables
template = \
u"""
---|------|------------------------|----|-------|-------|------|-----|----|---------|--------------|
 No Hit    Short Desc               Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
"""

def hhparse(hhresult_file, verbose):
	'''Convert HHpred's text-based output table in to an actual matrix'''
	pattern = StringIO(template).readlines()[1]
	colBreaks = [i for i, ch in enumerate(pattern) if ch == '|']
	widths = [j-i for i, j in zip( ([0]+colBreaks)[:-1], colBreaks ) ]

	hhtable = pd.read_fwf(hhresult_file, skiprows=8, nrows=10, header=0, widths = widths)
	if verbose is True:
		print(hhtable)

	top_hit = str(hhtable.loc[0,'Hit'])[0:4]
	top_hit_full = hhtable.loc[0,'Hit']
	top_prob = hhtable.loc[0,'Prob']
	top_eval = hhtable.loc[0,'E-value']
	top_pval = hhtable.loc[0,'P-value']
	top_score = hhtable.loc[0,'Score']

	if verbose is True:
		print("Your best hit: (PDB ID | Probability | E-Value | P-Value | Score)")
        	print("\t" +  str(top_hit) + "\t" + str(top_prob) + "\t" + str(top_eval) + "\t" + str(top_pval) + "\t" + str(top_score) )
	
	return top_hit, top_hit_full, top_prob, top_eval, top_pval, top_score

def getFullDesc(hhresult_file,top_hit_full, verbose):
	with open(hhresult_file, 'r') as hrh:
		for line in hrh:
			if line.startswith('>' + top_hit_full):
				full_desc = line

	if verbose is True:
		print full_desc

	return full_desc

def main():

	try:
   		parser = argparse.ArgumentParser(description='This script converts HHpreds verbose output in to a full table for subsequent analysis.')
		parser.add_argument(
			'-i',
			'--infile',
			action='store',
			required=True,
			help='The HHpred output file to parse.')
		parser.add_argument(
			'-v',
			'--verbose',
			action='store_true',
			help='Print additional messages to screen. Default behaviour is false, only the result would be printed to screen for piping etc.')
		parser.add_argument(
			'-o',
			'--outfile',
			action='store',
			default='None',
			help='Output file name to store results in. If none provided, the default will be infile.tsv.')			
					
		args = parser.parse_args()
		
	except:
		print "An exception occured with argument parsing. Check your provided options."
		traceback.print_exc()


	verbose = args.verbose
	hhresult_file = args.infile
	
	indir = os.path.dirname(hhresult_file)

	split = os.path.splitext(args.infile)
	basename = os.path.basename(split[0])

	if args.outfile is 'None':
		outfile = indir+ '/' + basename + '.tsv'
	else:
		outfile = args.outfile

# Main code begins:
	
	top_hit, top_hit_full, top_prob, top_eval, top_pval, top_score = hhparse(hhresult_file, verbose)
	full_desc = getFullDesc(hhresult_file, top_hit_full, verbose)

	
	with open(outfile, 'w') as ofh:
		ofh.write(basename + "\t" + str(top_hit) + "\t" + str(top_hit_full) + "\t" + str(top_prob) + "\t" + str(top_eval) + "\t" + str(top_pval) + "\t" + str(top_score) + "\t" + full_desc )

	print(basename + "\t" + str(top_hit) + "\t" + str(top_hit_full) + "\t" + str(top_prob) + "\t" + str(top_eval) + "\t" + str(top_pval) + "\t" + str(top_score) + "\t" + full_desc)



if __name__ == '__main__':
	main()
