# -*- coding: utf-8 -*-

"""
This script takes the .hhr files output by HHSuite and
turns the quite verbose file in to a fully tabulated
version with all the fields separated one, one line per
file. Thus, the file can be viewed simply in Excel etc.

It requires the non-standard pandas module.
"""

# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# The license for HHSuite remains with the authors. Please consult
# their licensce agreements before using this software, and ensure
# they are cited for any use. Run the script with -b to print rel-
# evant references.


import os
import sys
import traceback
import warnings

__author__ = "Joe R. J. Healey"
__version__ = "1.0.0"
__title__ = "tabulateHHpred"
__license__ = "GPLv3"
__author_email__ = "J.R.J.Healey@warwick.ac.uk"

template = \
u"""
---|------|------------------------|----|-------|-------|------|-----|----|---------|--------------|
 No Hit    Short Desc               Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
"""


def hhparse(hhresult_file, verbose):
	'''Convert HHpred's text-based output table in to a pandas dataframe'''
	from io import StringIO
	import pandas as pd

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
				full_desc = line.rstrip('\n')

	if verbose is True:
		print full_desc

	return full_desc

def getDOI(top_hit):
	"""Query the PDB REST API to get an associated DOI/Publication"""
	import requests

	query = requests.get("https://www.ebi.ac.uk/pdbe/api/pdb/entry/publications/" + str(top_hit))
	qjson = query.json()
	doi = qjson[top_hit][0]['doi']

	if not doi:
		doi = "No DOI found."
	return doi

def displayRefs():
    """Display relevant references"""
    print('''
    - The following publications pertain to HHSuite:
	- Alva V., Nam SZ., Söding J., Lupas AN. (2016)
	  The MPI bioinformatics Toolkit as an integrative platform for advanced protein sequence and structure analysis.
	  Nucleic Acids Res. pii: gkw348. PMID: 27131380
	
	- Remmert M., Biegert A., Hauser A., Söding J. (2011) HHblits: Lightning-fast iterative protein sequence searching by HMM-HMM alignment.
	  Nat Methods. 9(2):173-5. doi: 10.1038/nmeth.1818. PMID: 22198341
    
	- Söding J., Biegert A., Lupas AN. (2005) The HHpred interactive server for protein homology detection and structure prediction.
	  Nucleic Acids Res 33 (Web Server issue), W244-W248. PMID: 15980461
	
	- Söding J. (2005) Protein homology detection by HMM-HMM comparison.
	  Bioinformatics 21: 951-960. PMID: 15531603
	
	- Hildebrand A., Remmert A., Biegert A., Söding J. (2009) Fast and accurate automatic structure prediction with HHpred.
	  Proteins 77(Suppl 9):128-32. doi: 10.1002/prot.22499. PMID: 19626712
	
	- Meier A., Söding J. (2015) Automatic Prediction of Protein 3D Structures by Probabilistic Multi-template Homology Modeling.
	  PLoS Comput Biol. 11(10):e1004343. doi: 10.1371/journal.pcbi.1004343. PMID: 26496371
	''')
    sys.exit(1)

def parseArgs():
	"""Parse commandline arguments"""

	import argparse
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
		parser.add_argument(
			'-b',
			'--bibliography',
			action='store_true',
			help='Display references.')
							
		args = parser.parse_args()
		
	except:
		print "An exception occured with argument parsing. Check your provided options."
		traceback.print_exc()

	return parser.parse_args()

def main():
	"""Make function calls and handle arguments for tabulateHHpred.py"""
	args = parseArgs()

        if args.bibliography is True:
                displayRefs()

	verbose = args.verbose
	hhresult_file = args.infile
	
	indir = os.path.dirname(os.path.abspath(hhresult_file))

	split = os.path.splitext(args.infile)
	basename = os.path.basename(split[0])
	
	if args.outfile is 'None':
		outfile = indir + '/' + basename + '.tsv'
	else:
		outfile = args.outfile

# Main code begins:
	
	top_hit, top_hit_full, top_prob, top_eval, top_pval, top_score = hhparse(hhresult_file, verbose)
	full_desc = getFullDesc(hhresult_file, top_hit_full, verbose)
	doi = getDOI(top_hit)
	
	with open(outfile, 'w') as ofh:
		ofh.write(basename + "\t" + str(top_hit) + "\t" + str(top_hit_full) + "\t" + str(top_prob) + "\t" + str(top_eval) + "\t" + str(top_pval) + "\t" + str(top_score) + "\t" + full_desc + '\t' + doi)

	print(basename + "\t" + str(top_hit) + "\t" + str(top_hit_full) + "\t" + str(top_prob) + "\t" + str(top_eval) + "\t" + str(top_pval) + "\t" + str(top_score) + "\t" + full_desc + '\t' + doi)



if __name__ == '__main__':
	main()
