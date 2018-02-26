#!/home/wms_joe/bin/miniconda2/bin/python
"""
Match 2 PDB protein structures using UCSF Chimera on the commandline,
via Pychimera.
"""

import os
import subprocess
import sys
import argparse
import fnmatch

# If the script is invoked with a normal python interpreter,
# the environment must be patched. Not necessary if called
# with pychimera

# Ensure the shell environment variable $CHIMERADIR is set and
# points to the root directory of the Chimera installation.
# Should match the output of: $ chimera --root


if not sys.argv[0].endswith("pychimera"):
    import pychimera
    pychimera.patch_environ()
    pychimera.enable_chimera()

import chimera
from chimera import openModels, Molecule
from chimera import runCommand as rc
from MatchMaker import (match, CP_BEST,	GAP_OPEN,GAP_EXTEND, defaults, MATRIX, ITER_CUTOFF)

##### Main code begins #####


def main():
	try:
   		parser = argparse.ArgumentParser(description='Structurally match 1 or more query protein structures to a specified PDB entry via Pychimera and UCSF Chimera.')
		parser.add_argument(
			'-p',
			'--pdb',
			action='store',
			help='The 4 character PDB ID for the reference structure to match against.')
		parser.add_argument(
			'-s',
			'--subjects',
			nargs='+',
			help='A folder of PDB files to match to the reference structure.')
		parser.add_argument(
			'-o',
			'--outdir',
			action='store',
			help='Output directory for RMSD and session files.')
		args = parser.parse_args()

	except:
		print("An exception occured with argument parsing. Check your provided options.")
		sys.exit(1)

#TODO: add a switch for a reference PDB provided directly in case of a custom PDB.
	# Get reference structure from PDB
	print("\n")
	print("Beginning Chimera:")
	print("------------------")
	chimera.openModels.open(args.pdb,type="PDB")
	print("Opened reference structure: " + args.pdb)

	# Open model structures
	dirname = os.path.dirname(os.path.abspath(args.subjects[0]))
	for model in args.subjects:
		chimera.openModels.open(model, type="PDB")
		print("Successfully opened model: " + model)

	print("\n")
	print("Calculating RMSD (this may take a while):")
	print("-----------------------------------------")
	all_models = chimera.openModels.list(modelTypes=[chimera.Molecule])
	ref = all_models[0]
	sims = all_models[1:]

	for atoms1, atoms2, rmsd, fullRmsd in match(CP_BEST, [ref, sims], defaults[MATRIX], "nw", defaults[GAP_OPEN], defaults[GAP_EXTEND]):
		ref_mol = atoms2[0].molecule
		sim_mol = atoms1[0].molecule

		print(ref_mol.name + "\t" + sim_mol.name + "\t" + str(rmsd))

	with open("{}_RMSDs.tsv".format(os.path.join(dirname,os.path.basename(dirname))), 'w') as rmsd_fh:
		rmsd_fh.write(ref_mol.name + '\t' + sim_mol.name + '\t' + str(rmsd))

 	rc('save {0}_session.py'.format(os.path.join(dirname,os.path.basename(dirname))))
	print('Session saved.')
	chimera.closeSession()
	rc('stop now')
	print('Chimera exited. Results are stored in:' + dirname)
	sys.exit(0)

if __name__ == '__main__' :
	main()
