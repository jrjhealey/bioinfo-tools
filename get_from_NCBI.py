import argparse, subprocess, sys, fnmatch
from ftplib import FTP

#initiate parser
parser = argparse.ArgumentParser(description="This program will download bacterial sequencing data from NCBI's FTP server")

#optional arguments
parser.add_argument("-b", "--bacterium", help="name of bacterium (eg. Klebsiella_pneumoniae)", type=str) 
parser.add_argument("-d", "--database", default="refseq", choices = ["refseq", "genbank"], type=str)
parser.add_argument("-t", "--type", default="fasta", help="file type", choices = ["genbank", "fasta", "gff", "feature_table"], type=str)
parser.add_argument("-g", "--genomes", help="Will list all available genomes in specified database", choices=["genbank", "refseq"], type=str)
parser.add_argument("-o", "--outputpath", default=".", help="output directory", type=str)

#initialize FTP server
ftp_site = 'ftp.ncbi.nlm.nih.gov'
ftp = FTP(ftp_site)

#print help if no arguments given
if len(sys.argv) == 1:
	parser.print_help()
	sys.exit()
args = parser.parse_args()

#set default variables
NCBI_database = "refseq"
file_type = "fasta"
	
#change default variable based on command arguments (eg. -t genbank)
if args.database == "refseq":
	NCBI_database = "refseq"
elif args.database == "genbank":
	NCBI_database = "genbank"
	
if args.type == "fasta":
	file_type = "--exclude='*cds_from*' --exclude='*rna_from*' --include='*genomic.fna.gz' --exclude='*'"
elif args.type == "genbank":
	file_type = "--include='*genomic.gbff.gz' --exclude='*'"
elif args.type == "gff":
	file_type = "--include='*genomic.gff.gz' --exclude='*'"
elif args.type == "feature_table":
	file_type = "--include='*feature_table.txt.gz' --exclude='*'"
	
#Have a way to list all available genomes
if args.genomes:
	ftp.login()
	ftp.cwd('genomes/%s/bacteria' % NCBI_database)
	dirs = ftp.nlst()
	for genome in dirs:
		print(genome)
	sys.exit()
	
#Check that bacterium is in NCBI genbank/refseq
if args.bacterium:
	ftp.login()
	ftp.cwd('genomes/%s/bacteria' % NCBI_database)
	dirs = ftp.nlst()
	pattern = args.bacterium
	genome_match = fnmatch.filter(dirs, pattern)
	if len(genome_match) == 0:
		print("Bacterium %s not found in %s -----> Exiting program" % (args.bacterium, NCBI_database))
		sys.exit(1)
else:
	print("")
	print("Please provide bacterium name '-b'")
	print("Exiting program")
	print("")
	sys.exit(1)

#The meat of the program	
try:
	subprocess.call("rsync -Lrtv %s rsync://ftp.ncbi.nlm.nih.gov/genomes/%s/bacteria/%s/latest_assembly_versions/*/ %s" % (file_type, NCBI_database, args.bacterium, args.outputpath), shell=True)
except:
	raise Exception("Failed to download files")
