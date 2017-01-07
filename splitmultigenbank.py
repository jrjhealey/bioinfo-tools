# Usage:
# python splitmultigenbank.py multigenbank.gbk

# Note, the script will just dump out all the GBKs in to the same
# folder, named according to their record IDs.

from Bio import SeqIO
import sys

for rec in SeqIO.parse(sys.argv[1], "genbank"):
   SeqIO.write([rec], open(rec.id + ".gbk", "w"), "genbank")
