from ete3 import NCBITaxa
import ncbi_genome_download as ngd
import sys

from BioRanges.lightweight import Range, SeqRange
  # See also: bx-python, Quicksect and GenomicRanges



taxon_name = sys.argv[1]

ncbi = NCBITaxa()
ncbi.update_taxonomy_database()
ebact = ncbi.get_descendant_taxa(taxon_name)

for i in ebact:
	print("Fetching TAXID: " + str(i) + " ----> " +  ncbi.translate_to_names(i))
#	ngd.download(section='refseq', taxid=i)
