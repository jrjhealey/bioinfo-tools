# Usage:
# python splitmultigenbank.py multigenbank.gbk

# Note, the script will just dump out all the GBKs in to the same
# folder, named according to their record IDs.

from Bio import SeqIO
import sys
import argparse
import traceback
import warnings


def main():

# Parse command line arguments
        try:
                parser = argparse.ArgumentParser(
                        description='This script separates out multi-genbank files into individual ones')
                parser.add_argument(
                        '-g',
                        '--genbank',
                        action='store',
                        required=True,
                        help='The genbank file you wish to split.')

                args = parser.parse_args()

        except:
                print "An exception occured with argument parsing. Check your provided options."
                traceback.print_exc()

	for rec in SeqIO.parse(genbank, "genbank"):
  		SeqIO.write([rec], open(rec.id + ".gbk", "w"), "genbank")


if __name__ == '__main__' :
	main()
