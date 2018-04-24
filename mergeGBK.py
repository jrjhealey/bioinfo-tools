
# Merge genbanks into a single entry.
# E.g. for removal of multiple contigs


from Bio import SeqIO
import sys, os

file = sys.argv[1]
filename = os.path.splitext(sys.argv[1])[0]
outfile = sys.argv[2]

recs = list(SeqIO.parse(file, 'genbank'))

with open(outfile, 'w') as ofh:
	for i in recs:
		i.description = filename
		i.annotations = filename
		i.name = filename
		i.accession = filename
		SeqIO.write(i, ofh, 'genbank')

