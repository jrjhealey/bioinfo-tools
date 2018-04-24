

from Bio import SeqIO
import sys
import argparse
import traceback
import warnings


dinucleotides = ['AA','AT','AC','AG',
               'TT','TA','TC','TG',
               'CC','CA','CT','CG',
               'GG','GA','GT','GC']

ratio = ['AA/TT','AC/GT','AG/CT',
         'TT/AA','TC/GA','TG/CA',
         'CC/GG','CA/TG','CT/AG',
         'GG/CC','GA/TC','GT/AC']


def parseArgs():
	"""Parse commandline arguments"""

	import argparse
	try:
   		parser = argparse.ArgumentParser(description='Count dinucleotide frequency in a multifasta and write the output to a CSV.')
		parser.add_argument('--infile',
				    action='store',
				    help='The HHpred output file to parse.')
		parser.add_argument('--ratiofile',
				    action='store',
				    help='Output file to store ratios in.')
		parser.add_argument('--countsfile',
				    action='store',
				    help='Output file to store counts in.')	
	except:
		print "An exception occured with argument parsing. Check your provided options."
		traceback.print_exc()

	return parser.parse_args()



def Ratio(c):
	"""This function calculates the ratio between each dinucleotide on the sense and antisense."""
	di = dict() #A dictionary for the ratio.

	for ra in ratio:
		r = ra.split('/') #Splits the header for 2 dinucleotides.(r is a list)
		first = r[0] #The first dinucleotide.
		second = r[1] #The second dinucleotide.
		if(c[second] == 0): #In case of division at 0.
			di[ra] = None
	else:
		di[ra] = ((float)(c[first])/c[second]) #The ratio beetween the first dinucleotide and the second.

	return di


def Counting(seq):
        """This function gets a cDNA and returns the frequency of each dinucleotide"""
        counting = {k: 0 for k in dinucleotides}

        for i in range(len(seq)-2):
                if seq[i:i+2] in counting:
                        counting[seq[i:i+2]] += 1

        return counting



def main():
	import csv

	args = parseArgs()

	with open(args.infile, 'r') as ifh, open(args.ratiofile, 'wb') as rfh, open(args.countsfile, 'wb') as cfh:
		for record in SeqIO.parse(ifh, 'fasta'):
			identifier = record.id
			sequence = record.seq
			sequence = sequence.upper()

			count = Counting(sequence)

			c = csv.DictWriter(cfh, count.keys())
			cfh.write(identifier + ',')
			c.writerow(count)

			ratio = Ratio(count)	
			r = csv.DictWriter(rfh, ratio.keys())
			rfh.write(identifier + ',')
			r.writerow(ratio)




if __name__ == '__main__':
	main()
