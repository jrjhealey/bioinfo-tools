# Usage:
# python extractNAfeatures.py genbank.gbk outfile.ffn

from Bio import SeqIO
import sys

with open(sys.argv[2], 'w') as nfh:
        for rec in SeqIO.parse(sys.argv[1], "genbank"):
                if rec.features:
                        for feature in rec.features:
                                if feature.type == "CDS":
                                        nfh.write(">%s from %s\n%s\n" % (
                                        feature.qualifiers['gene'][0],
                                        rec.name,
                                        feature.location.extract(rec).seq))
