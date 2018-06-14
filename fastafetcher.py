
# Extract fasta files by their descriptors stored in a separate file.
# Requires biopython

from Bio import SeqIO
import sys, six
import argparse

def getKeys(args):
	"""Turns the input key file into a list. May be memory intensive."""

	with open(args.keyfile, "r") as kfh:
		keys = []
    		for line in kfh:
			line = line.rstrip('\n')
			line = line.lstrip('>')
			keys.append(line)

	return keys

def main():
	"""Takes a string or list of strings in a text file (one per line) and retreives them and their sequences from a provided multifasta."""
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
			'--keyfile',
			action='store',
			help='A file of header strings to search the multifasta for. Must be exact. Must be one per line.')
		parser.add_argument(
			'-s',
			'--string',
			action='store',
			help='Provide a string to look for directly, instead of a file.')
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
		print('An exception occured with argument parsing. Check your provided options.')
		sys.exit(1)
# Main code:
	# Call getKeys() to create the list of keys from the provided file:
	
	if not (args.keyfile or args.string):
		print('No key source provided. Exiting.')
		sys.exit(1)
	elif not args.keyfile:
		keys = args.string
	else:
		keys = getKeys(args)
	
	if args.verbose is not False:
		if args.invert is False:
			print('Fetching the following keys: ' + args.keyfile)
			if isinstance(keys, six.string_types):
				print(keys)
			elif isinstance(keys, (list,)):
				for key in keys:
					print(key)
		else:
			print('Ignoring the following keys, and retreiving everything else from: ' + args.fasta)
			for key in keys:
				print(key)

	# Parse in the multifasta and assign an iterable variable:
        seqIter = SeqIO.parse(args.fasta, 'fasta')

	# For each sequence in the multifasta, check if it's in the keys[] tuple. If so, print it out:
	for seq in seqIter:
		if args.invert is False:
			if seq.id in keys:
				print(seq.format("fasta"))
			if args.outfile is not None:
				SeqIO.write(seq, outfile, "fasta")
		else:
	# If the --invert/-i flag is given, print all fastas NOT listed in the keyfile
			if seq.id not in keys:
				print(seq.format("fasta"))
			if args.outfile is not None:
				SeqIO.write(seq, outfile, "fasta")

if __name__ == "__main__":
	main()

