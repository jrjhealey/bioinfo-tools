#!/home/wms_joe/bin/miniconda2/bin/python
"""
This script pulls in homologs of proteins from PDB
as determined by HHSuite. It then employs UCSF Chimera to
structurally match them and get an indication of how well
they score (RMSD) in order to pick the best simulation.
"""
import os
import subprocess
import sys
import argparse
import traceback
import warnings
import fnmatch
import pandas as pd
from io import StringIO

import pychimera

import chimera
from chimera import openModels, Molecule
from chimera import runCommand as rc
from MatchMaker import (match,
 						CP_BEST,
						GAP_OPEN,
						GAP_EXTEND,
						defaults,
						MATRIX,
						ITER_CUTOFF)

# HHSuite output parser/template

template = \
u"""
---|------|------------------------|----|-------|-------|------|-----|----|---------|--------------|
 No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
"""

def hhparse(hhresult_file):
    pattern = StringIO(template).readlines()[1]
    colBreaks = [i for i, ch in enumerate(pattern) if ch == '|']
    widths = [j-i for i, j in zip( ([0]+colBreaks)[:-1], colBreaks ) ]

    hhtable = pd.read_fwf(hhresult_file, skiprows=8, nrows=10, header=0, widths = widths)
    print(hhtable.dtypes)
    print(hhtable)

    top_hit = str(hhtable.loc[0,'Hit'])[0:4]
    top_prob = hhtable.loc[0,'Prob']
    top_eval = hhtable.loc[0,'E-value']
    top_pval = hhtable.loc[0,'P-value']
    top_score = hhtable.loc[0,'Score']

    return top_hit, top_prob, top_eval, top_pval, top_score

# If running using python interpreter and not pychimera:
# os.environ['CHIMERADIR'] = '/home/wms_joe/Applications/CHIMERA1.11'
# CHIMERADIR should point to the application root directory.
# This can be found with:   `chimera --root`


# Only needed if running via normal python interpreter, not pychimera
# pychimera.patch_environ()
# pychimera.enable_chimera()


##### Main code begins #####



# Step 1, parse the output of HHpred to get the nearest homolog.

def main():
	def_cpus=4
	def_dir=os.getcwd()
	try:
   		parser = argparse.ArgumentParser(description='This script compares protein structural homologs as determined with HHpred, to PDB models using UCSF CHIMERA, to gather metrics of structural similarity.')

		parser.add_argument(
			'-d',
			'--database',
			action='store',
			default=None,
			help='You can specify a different HHpred database (filepath) to use if you offer this parameter. Otherwise it defaults to PDB.')

		parser.add_argument(
			'-f',
			'--fasta',
			action='store',
			help='The fasta amino acid sequence that corresponds to the simulated structure and the sequence you wish to query.')

		parser.add_argument(
			'-s',
			'--simulations',
			action='store',
			help='The directory where the protein structure simulations are stored. They must be in PDB format and be named logically.')

		parser.add_argument(
			'-r',
			'--rmsd',
			action='store',
			default=None,
			help='The filename to store the RMSD values for each matching.')

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
		split = os.path.splitext(args.fasta)
		basename = os.path.basename(split[0])
	else:
		print('No input fasta was provided. Check your input parameters')
		sys.exit(1)

	# Run the HHsearch
	if args.fasta is not None and args.database is not None:
		print(sys.path)
		print("\n")
		print("Running " + args.fasta + " in HHsearch, against " + args.database + " on " + str(args.threads) + " threads.")
		hhresult_file = './{0}.hhr'.format(basename)
		search_cmd = 'hhsearch -dbstrlen 50 -B 1 -b 1 -p 60 -Z 1 -E 1E-03 -nocons -nopred -nodssp -cpu {0} -i {1} -d {2}'.format(args.threads,args.fasta, args.database)
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

        top_hit, top_prob, top_eval, top_pval, top_score = hhparse(hhresult_file)
        print("Your best hit: (PDB ID | Probability | E-Value | P-Value | Score)")
        print("\t" +  str(top_hit) + "\t" + str(top_prob) + "\t" + str(top_eval) + "\t" + str(top_pval) + "\t" + str(top_score) )

	print("Full results are found in: " + hhresult_file)

	# Acquire the protein simulations

	model_list = []
	for root, dirnames, filenames in os.walk(args.simulations):
		for filename in sorted(fnmatch.filter(filenames, '*.pdb')):
			if filename.startswith(basename):
				model_list.append(os.path.join(root, filename))

	if not model_list:
		print("No protein structures were found matching: " + basename)
		sys.exit(1)
	else:
		print("\n")
		print("Found the following models:")
		print("---------------------------")
		for model_path in model_list:
			locus_dir = os.path.dirname(os.path.abspath(model_path))
			print("Found: " + os.path.basename(model_path) + " in " + locus_dir)


	# Get reference structure from PDB
	print("\n")
	print("Beginning Chimera:")
	print("---------------------------")
	chimera.openModels.open(top_hit,type="PDB")
	print("Opened reference structure: " + top_hit)

	# Open model structures
	for model_path in model_list:
		chimera.openModels.open(model_path,type="PDB")
		print("Opened: " + os.path.basename(model_path))

 	if args.rmsd is not None:
		rmsd_file = args.rmsd
	else:
		joint_path = os.path.join(locus_dir,basename)
 		rmsd_file = '{0}_RMSDs.tsv'.format(joint_path)

	print("\n")
	print("RMSD values are:")
	print("---------------------------")
 	with open(rmsd_file, "w") as rmsd_output_file:
		all_models = chimera.openModels.list(modelTypes=[chimera.Molecule])
		ref = all_models[0]
		sims = all_models[1:]

		for atoms1, atoms2, rmsd, fullRmsd in match(CP_BEST,[ref, sims],defaults[MATRIX],
													"nw",defaults[GAP_OPEN],defaults[GAP_EXTEND]):
			ref_mol = atoms2[0].molecule
			sim_mol = atoms1[0].molecule

			print(ref_mol.name + "\t" + sim_mol.name + "\t" + str(rmsd))
			rmsd_output_file.write(ref_mol.name + "\t" \
                                    + sim_mol.name + "\t" \
                                    + str(rmsd) + "\t" \
                                    + str(top_prob) + "\t" \
                                    + str(top_eval) + "\t" \
                                    + str(top_pval) + "\t" \
                                    + str(top_score) + "\n")

 	rc('save {0}_session.py'.format(joint_path))
	chimera.closeSession()
	print("Chimera exited. All done. All results are in:")
	print(locus_dir)
	rc('stop now')
	sys.exit(0)

if __name__ == '__main__' :
	main()
