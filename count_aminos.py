import collections
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import sys

for rec in SeqIO.parse(sys.argv[1], "fasta"):
    x = ProteinAnalysis(str(rec.seq))
#    if sys.argv[2] is "sort":
    for key, val in sorted(x.iteritems(), key=lambda (k,v): (v,k)):
        print "%s %s:%s" % (key, value)
#    else:
#        print rec.id, x.count_amino_acids()
