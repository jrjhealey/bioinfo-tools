#!/usr/bin/python

# Script to parse orthoMCL/orthAgogue gene clusters and extract just those that are core

import re
import sys
import argparse
import traceback
import warnings

##################

# Get a tuple of strings to iterate over 
def getKeys(nameFile):
	with open(nameFile, "r") as namehandle:
		names = []
    		for line in namehandle:
	        	strip = line.rstrip('\n')
			names.append(strip)
        
	return names

#################
def main():
###################################################################################################
# Parse command line arguments
	try:
   		parser = argparse.ArgumentParser(description='This script slices entries such as genes or operons out of a genbank, subsetting them as their own file.')
		parser.add_argument(
			'-i',
			'--infile',
			action='store',
			help='The gene cluster file to parse')
		parser.add_argument(
			'-n',
			'--names',
			action='store',
			help='A file of compliant "keys" which correspond to the genomes from which the clusters are derived.')
		args = parser.parse_args()
	except:
		print "An exception occured with argument parsing. Check your provided options."
		traceback.print_exc()

	nameFile = args.names

# Main code:
	keys = getKeys(nameFile)
	print(keys)
	matchedLines = []
	with open(args.infile, "r") as clusterFile:
		matchedLines = [line for line in clusterFile if all([line.count(key) == 1 for key in keys])]

	for each in matchedLines:
		print each.rstrip('\n')

if __name__ == "__main__":
	main()

