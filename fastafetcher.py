
# Extract fasta files by their descriptors stored in a separate file.

from Bio import SeqIO
import sys
import traceback
import warnings
import argparse

def getKeys(keyFile):
	"""Turns the input key file into a tuple. May be memory intensive."""

	with open(keyFile, "r") as kfh:
		keys = []
    		for line in kfh:
			line = line.rstrip('\n')
			line = line.lstrip('>')
			keys.append(line)

	return keys

def main():
	"""Takes a list of strings in a text file (one per line) and retreives them and their sequences from a provided multifasta."""
	# Parse arguments from the commandline:
	try:
   		parser = argparse.ArgumentParser(description='Retrieve one or more fastas from a given multifasta.')
		parser.add_argument(
			'-i',
			'--infile',
			action='store',
			help='The multifasta to search.')
		parser.add_argument(
			'-k',
			'--keys',
			action='store',
			help='A file of header strings to search the multifasta for. Must be exact. Must be one per line.')
		parser.add_argument(
			'-o',
			'--outfile',
			action='store',
			default=None,
			help='Output file to store the new fasta sequences in. Just prints to screen by default.')
		parser.add_argument(
			'-v',
			'--verbose',
			action='store_true',
			help='Set whether to print the key list out before the fasta sequences. Useful for debugging.')

		args = parser.parse_args()
	except:
		print "An exception occured with argument parsing. Check your provided options."
		traceback.print_exc()

	# Rename args to enable them to be provided to the getKeys function:
	keyFile = args.keys
	inFile = args.infile

# Main code:
	# Call getKeys() to create the tuple of keys from the provided file:
	keys = getKeys(keyFile)

	if args.verbose is not False:
		print('Fetching the following keys from: ' + inFile)
		for key in keys:
			print(key)

	# Parse in the multifasta and assign an iterable variable:
        seqIter = SeqIO.parse(inFile, 'fasta')

	# For each sequence in the multifasta, check if it's in the keys[] tuple. If so, print it out:
	for seq in seqIter:
		if seq.id in keys:
			print(seq.format("fasta"))

#	SeqIO.write((seq for seq in seqIter if seq.id in keys), sys.stdout, "fasta")

if __name__ == "__main__":
	main()

