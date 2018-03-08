from ete3 import NCBITaxa
import ncbi-genome-download as ngd
import sys

taxon_name = sys.argv[1]

ncbi = NCBITaxa()
ncbi.update_taxonomy_database()
ebact = ncbi.get_descendant_taxa(taxon_name)

for i in ebact:
	ngd.download(section='refseq', group='bacteria', taxid=i)
