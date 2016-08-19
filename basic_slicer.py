#!/usr/bin/python

# This script is designed to take a genbank file and 'slice out'/'subset'
# regions (genes/operons etc.) and produce a separate file.

# Based upon the tutorial at http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc44

# Next dev task, implement reading column of co-ords from a blast tabular result?
# Set up and handle arguments:
from Bio import SeqIO
import os, getopt, sys


def main(argv):
	record = ''
	start = ''
	end = ''
	try:
		opts, args = getopt.getopt(argv, 'hi:o:s:e:', [
													   'help',
													   'input=',
													   'outfile=',
													   'start=',
													   'end='
													   ]
								  )
		if not opts:
	 		print "No options supplied. Aborting."
	 		usage()
	 		sys.exit(2)
	except getopt.GetoptError:
		print "Some issue with commandline args.\n"
		usage()
		sys.exit(2)

	for opt, arg in opts:
		if opt in ("-h", "--help"):
			usage()
			sys.exit(2)
		elif opt in ("-i", "--input"):
			filename = arg
			record = SeqIO.read(arg, "genbank")
		elif opt in ("-o", "--outfile"):
			outfile = arg
		elif opt in ("-s", "--start"):
			start = int(arg)
		elif opt in ("-e", "--end"):
			end = int(arg)
	print("Slicing " + filename + " from " + str(start) + " to " + str(end))
	sub_record = record[start:end]
	SeqIO.write(sub_record, outfile, "genbank")


def usage():
	print(
"""
This script 'slices' entries such as genes or operons out of a genbank,
subsetting them as their own file.\n

Usage:
python slice_genbank.py -h|--help -i|--input <genbank> -o|--output <genbank> -s|--start <int> -e|--end <int>"

Options:

	-h|--help		Displays this usage message. No options will also do this.
	-i|--input		The genbank file you which to subset a record from.
	-o|--outfile	The file name you wish to give to the new sliced genbank.
	-s|--start		An integer base index to slice the record from.
	-e|--end		An integer base index to slice the record to.
"""
		  )

if __name__ == "__main__":
   main(sys.argv[1:])
