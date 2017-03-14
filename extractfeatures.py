# Usage:
# python extractfeatures.py infile.gbk outfile.faa

from Bio import SeqIO
import sys

with open(sys.argv[2], 'w') as ofh:
        for seq_record in SeqIO.parse(sys.argv[1], 'genbank'):
                for seq_feature in seq_record.features:
                        if seq_feature.type=="CDS":
                                assert len(seq_feature.qualifiers['translation'])==1
                                ofh.write(">%s from %s\n%s\n" % (
                                seq_feature.qualifiers['gene'][0],
                                seq_record.name,
                                seq_feature.qualifiers['translation'][0]))
