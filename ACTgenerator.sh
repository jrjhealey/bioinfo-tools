#!/bin/bash

# A small script to generate artemis comparison files (nucleic acid comparison)
# since all the webservers are apparently defunct!

# Script requires blastn (NOT LEGACY BLAST)  and makeblastdb in path - check for existence:

command -v makeblastdb >/dev/null 2>&1 || { echo >&2 "makeblastdb doesn't appear to be installed. Aborting."; exit 1; }
command -v blastn >/dev/null 2>&1 || { echo >&2 "BLAST+ doesn't appear to be installed. Aborting."; exit 1; }
# Capture inputs

usage()
{
cat << EOF
usage: $0 options

This script generates the necessary BLAST comparison file
for use with the Artemis Comparison tool, when comparing 2 genomes.

OPTIONS:
   -h | --help     	Show this message
   -r | --ref		Reference fasta sequence
   -q | --query		Query fasta sequence
   -d | --database	Name of the BLAST database comparison file.
   -o | --outdir	Directory to output all the files to
   -t | --tidy		Make the script tidy up after itself (T/F).
			Default F.

EOF
}

# Tolerate long arguments
for arg in "$@"; do
  shift
  case "$arg" in
  		"--help")
		set -- "$@" "-h"
		;;
	    	"--reference")
    		set -- "$@" "-r"
		;;
   		"--query")
    		set -- "$@" "-q"
    		;;
    		"--database")
       		set -- "$@" "-d"
       		;;
		"--outdir")
		set -- "$@" "-o"
		;;
		"--tidy")
		set -- "$@" "-t"
		;;
	    	*)
		set -- "$@" "$arg"
  esac
done
# getopts assigns the arguments to variables
while getopts "hr:q:d:o:t:" OPTION ;do
  case $OPTION in
		r)
		reference=$OPTARG
		;;
		q)
		query=$OPTARG
		;;
		d)
		database=$OPTARG
		;;
		o)
		outdir=$OPTARG
		;;
		t)
		tidy=$OPTARG
		;;
		h)
		usage
		exit
		;;
	esac
done

if [[ -z $reference ]]; then
    	usage
    	echo "No Reference sequence provided. Exiting." ; exit 1
fi

if [[ -z $query ]]; then
    	usage
    	echo "Query not supplied. Exiting." ; exit 1
fi

if [[ -z $database ]]; then
    	database=${reference%.*}_DB
	echo "No database name was specified. Using the reference sequence name and appending _DB."
fi

if [[ -z $outdir ]]; then
	outdir=$(pwd)
	echo "No output directory was specified. Defaulting to the current working directory:"
	pwd
fi

if [ $tidy == "True" ]; then
 tidy="T"
 else
  if [ $tidy == "False" ] || [ -z $tidy ]; then
  tidy="F"
  fi
fi



#####

# Step 1: Make a BLAST database of the reference sequence:

echo -en Executing command with arguments: '\n' Reference - $reference '\n' Database - $database '\n' Query sequence - $query '\n'

makeblastdb -in $reference -dbtype 'nucl' -title $database -out ${outdir%/}/${database} -parse_seqids
echo "Database created."

# Step 2: Perform the all-vs-all BLAST using the query sequence and reference database.
blastn -db $database -query $query -outfmt 6 -out ${outdir%/}/${query%.*}_vs_${reference%.*}.act

echo "All finished! The comparion file is called:"
echo "${query%.*}_vs_${reference%/*}.act"

if [ $tidy == "T" ]; then
	rm -v ${outdir%/}/${database}*
fi
