#!/bin/bash

# Script to align fastq reads to a given reference, and output the necessary files sorted etc.

# Script requires bwa - check for existence:

command -v bwa >/dev/null 2>&1 || { echo >&2 "BWA doesn't appear to be installed. Aborting."; exit 1; }

# Capture inputs

usage()
{
cat << EOF
usage: $0 options

This script aligns fastq reads to a given
reference, and outputs alignment maps etc.

OPTIONS:
   -h | --help     	   Show this message
   -r | --ref		   Reference fasta sequence
   -1 | --read1 	   First set of reads ("R1")
   -2 | --read2		   Second set of reads ("R2")

EOF
}

# Tolerate long arguments
for arg in "$@"; do
  shift
  case "$arg" in
  		"--help")
		set -- "$@" "-h"
		;;
    	"--ref")
    		set -- "$@" "-r"
		;;
   		"--read1")
    		set -- "$@" "-1"
    		;;
    	"--read2")
       		set -- "$@" "-2"
       		;;
    	*)
			set -- "$@" "$arg"
  esac
done
# getopts assigns the arguments to variables
while getopts "hr:1:2:" OPTION
do
	case $OPTION in
		r)
		reference=$OPTARG
		;;
		1)
		read1=$OPTARG
		;;
		2)
		read2=$OPTARG
		;;
		h)
		usage
		exit
		;;
	esac
done

if [[ -z $reference ]]
	then
    	usage
    	echo "No Reference sequence provided. Exiting." ; exit 1
fi
if [[ -z $read1 ]]
	then
    	usage
    	echo "Forward reads not supplied. Exiting." ; exit 1
fi
if [[ -z $read2 ]]
	then
    	usage
    	echo "Reverse reads not supplied. Exiting." ; exit 1
fi

reads=("$read1" "$read2")

# Index the reference sequence
bwa index "$reference"
echo "BWA FINISHED INDEXING THE REFERENCE"

# Align reads
alignedreads=()

for readelement in ${reads[@]} ; do
	bwa aln "$reference" "$readelement" > ${readelement%.*.*}.sai
	echo "BWA FINISHED ALIGNING $readelement "
	alignedreads+=("${readelement%.*.*}.sai")
done
echo "BWA FINISHED ALIGNING EACH READ."

bwa sampe "$reference" "${alignedreads[0]}" "${alignedreads[1]}" "$read1" "$read2" > "${reference%.*}"_aligned.sam
echo "SAM FILE CREATED"

samtools view -S "${reference%.*}"_aligned.sam -b -o ./"${reference%.*}"_aligned.bam
echo "SAM FILE CONVERTED TO BAM"

samtools sort "${reference%.*}"_aligned.bam "${reference%.*}"_aligned_sorted
echo "BAM FILE SORTED"

samtools index "${reference%.*}"_aligned_sorted.bam
echo "BAM FILE INDEXED"



