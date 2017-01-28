from Bio import SeqIO
import sys
import argparse
import traceback
import warnings
import os



def main():
###################################################################################################
# Parse command line arguments
        try:
                parser = argparse.ArgumentParser(
                        description='This script slices entries such as genes or operons out of a genbank, subsetting them as their own file.')
                parser.add_argument(
                        '-g',
                        '--genbank',
                        action='store',
                        required=True,
                        help='The genbank file to extract features from.')
                parser.add_argument(
                        '-o',
                        '--outfile',
                        action='store',
                        help='Fasta file to output to - the script will add extensions automatically for .ffn and .ffa.')
		parser.add_argument(
			'-m',
			'--method',
			action='store',
			default='both',
			help='Which feature format to output: "prot|nuc|both" (default "both").')

                args = parser.parse_args()

        except:
                print "An exception occured with argument parsing. Check your provided options."
                traceback.print_exc()

	genbank = args.genbank
	method = args.method

# Main code:
if args.method is "prot":
	with open(outfile + '.ffa', 'w') as pfh:
        	for seq_record in SeqIO.parse(args.genbank, 'genbank'):
                	for seq_feature in seq_record.features:
                        	if seq_feature.type=="CDS":
                                	assert len(seq_feature.qualifiers['translation'])==1
                                	ofh.write('>' +	seq_feature.qualifiers['gene'][0] + ' | ' + seq_record.name + '\n' + seq_feature.qualifiers['translation'][0] + '\n')

if args.method is "nucl":
	with open(sys.argv[2], 'w') as nfh:
        	for rec in SeqIO.parse(sys.argv[1], "genbank"):
                	if rec.features:
                        	for feature in rec.features:
                                	if feature.type == "CDS":
                                        	nfh.write(">%s from %s\n%s\n" % (
                                        	feature.qualifiers['gene'][0],
                                        	rec.name,
                                        	feature.location.extract(rec).seq))

def generateRecord(genbank)
	record = SeqIO.parse(genbank, "genbank")
	return record

def generateFasta(locus
	content = '>' + locus

if __name__ == "__main__":
        main()
