#!/usr/bin/python

"""
Script to get secondary structure proportions from a PDB
via chimera/pychimera

Usage:

pychimera SSCalc.py -i file.pdb [ -o path/to/outfile ]
"""
import os
import subprocess
import sys
import argparse
import traceback
import warnings
import fnmatch

import pychimera

import chimera
from chimera import openModels, Molecule
from chimera import runCommand as rc

# If running using python interpreter and not pychimera:
#if argv[0] == "pychimera" :
#	continue
#elif argv[0] == "python" :
#	os.environ['CHIMERADIR'] = '/home/wms_joe/Applications/CHIMERA1.11'
#	pychimera.patch_environ()
#	pychimera.enable_chimera()

# CHIMERADIR should point to the application root directory.
# This can be found with:   `chimera --root`
# Only needed if running via normal python interpreter, not pychimera

##### Main code begins #####

def main():
	try:
   		parser = argparse.ArgumentParser(description='Calculate PDB secondary structure proportions')

		parser.add_argument(
			'-i',
			'--infile',
			action='store',
			default=None,
			required=True,
			help='Input PDB to be analysed')

		parser.add_argument(
			'-o',
			'--outfile',
			action='store',
			default=None,
			help="Output filename. Default is the input file with the extension .ss appended in the current directory.")

		args = parser.parse_args()

	except:
		print("An exception occured with argument parsing. Check your provided options.")
		traceback.print_exc()

	#### Calculating ####
	if args.infile is not None:
		print('Calculating Secondary structure for ' + args.infile)
		basename = os.path.basename(os.path.splitext(args.infile)[0])
	else:
		print('No structure file was provided; exiting.')
		sys.exit(1)
	
	chimera.openModels.open(args.infile,type="PDB")
	if args.outfile is None:
		args.outfile = './{}.ss'.format(basename)

	print('Secondary Structure Proportions (AA Length | Helix % | Sheet % | Other %)')
	with open(args.outfile, "w") as outputFile:
		for mol in openModels.list(modelTypes=[Molecule]):
			helix_fract = len([r for r in mol.residues if r.isHelix]) / float(len(mol.residues))
			sheet_fract = len([r for r in mol.residues if r.isSheet]) / float(len(mol.residues))
			other_fract = (1 - (helix_fract + sheet_fract))
			helix_perc = helix_fract * 100
			sheet_perc = sheet_fract * 100
			other_perc = other_fract * 100
			print(args.infile + "\t" +
				 str(len(mol.residues)) + "\t" +
				 str(helix_perc) + "\t" +
				 str(sheet_perc) + "\t" +
				 str(other_perc))
			outputFile.write(args.infile + "\t" +
				 str(len(mol.residues)) + '\t' +
				 str(helix_perc) + "\t" +
				 str(sheet_perc) + "\t" +
				 str(other_perc) + "\n")

	chimera.closeSession()
	print("Chimera exited. All done. All results are in: " + args.outfile)
	rc('stop now')
	sys.exit(0)







if __name__ == '__main__' :
	main()
