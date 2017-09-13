# Order a multifasta by the order of sequence names in a phylip

def main():
	"""Order a fasta by the corresponding order of names in a phylip"""
	import argparse
	try:
   		parser = argparse.ArgumentParser(description='Output a multifasta ordered by an input file (phylip alignment).')
		parser.add_argument(
			'-f',
			'--fasta',
			action='store',
			required=True,
			help='The multifasta to re-rorder.')
		parser.add_argument(
			'-p',
			'--phylip',
			action='store',
			required=True,
			help='The phylip file which specifies the order of the sequences.')
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
		print('An exception occured with argument parsing. Check your provided options.')
		traceback.print_exc()

	from Bio import AlignIO, SeqIO

	align = AlignIO.read(args.phylip,'phylip')
	fasta_dict = SeqIO.to_dict(SeqIO.parse(args.fasta, 'fasta'))

	phynames = [rec.id for rec in align]

	for key in phynames:
        	print(fasta_dict[key].format('fasta'))
		if outfile is not None:
			SeqIO.write(fasta_dict[key], args.outfile, "fasta")


if __name__ == '__main__':
	main()
