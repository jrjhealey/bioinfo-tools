
# Extract fasta files by their descriptors stored in a separate file.
# Requires biopython


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
			'-f',
			'--fasta',
			action='store',
			required=True,
			help='The multifasta to search.')
		parser.add_argument(
			'-k',
			'--keys',
			action='store',
			required=True,
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
		parser.add_argument(
			'-i',
			'--invert',
			action='store_true',
			help='Invert the search, and retrieve all sequences NOT specified in the keyfile.')
		args = parser.parse_args()
	
	except:
		print "An exception occured with argument parsing. Check your provided options."
		traceback.print_exc()

	# Rename args to enable them to be provided to the getKeys function:
	keyFile = args.keys
	inFile = args.fasta
	outFile = args.outfile
# Main code:
	# Call getKeys() to create the tuple of keys from the provided file:
	keys = getKeys(keyFile)

	if args.verbose is not False:
		if args.invert is False:
			print('Fetching the following keys from: ' + inFile)
			for key in keys:
				print(key)
		else:
			print('Ignoring the following keys, and retreiving everything else from: ' + inFile)
			for key in keys:
				print(key)

	# Parse in the multifasta and assign an iterable variable:
        seqIter = SeqIO.parse(inFile, 'fasta')

	# For each sequence in the multifasta, check if it's in the keys[] tuple. If so, print it out:
	for seq in seqIter:
		if args.invert is False:
			if seq.id in keys:
				print(seq.format("fasta"))
			if args.outfile is not None:
				SeqIO.write(seq, outFile, "fasta")
		else:
	# If the --invert/-i flag is given, print all fastas NOT listed in the keyfile
			if seq.id not in keys:
				print(seq.format("fasta"))
			if args.outfile is not None:
				SeqIO.write(seq, outFile, "fasta")

if __name__ == "__main__":
	main()

