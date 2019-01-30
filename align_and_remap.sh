#!/bin/bash
set -eo pipefail

# Script to align fastq reads to a given reference, and output the necessary files sorted etc.

if [ -t 1 ] ; then
  ncols=$(tput colors)
  if [ -n "$ncols" ] && [ "$ncols" -ge 8 ] ; then
    bold="$(tput bold)"
    underline="$(tput smul)"
    rmunderline="$(tput rmul)"
    standout="$(tput smso)"
    black="$(tput setaf 0)"
    red="$(tput setaf 1)"
    green="$(tput setaf 2)"
    yellow="$(tput setaf 3)"
    blue="$(tput setaf 4)"
    magenta="$(tput setaf 5)"
    cyan="$(tput setaf 6)"
    white="$(tput setaf 7)"
    default="$(tput sgr0)"
  fi
fi

log(){
  # Logging function (prints to STDOUT in WHITE).
  echo -e >&1 "${white}${underline}INFO:${rmunderline} ${1:-$(</dev/stdin)}${default}"
}

err(){
  # Error function (prints to STDERR in RED).
  echo -e >&2 "${red}${underline}ERROR:${rmunderline} ${1:-$(</dev/stdin)}${default}"
}

warn(){
  # Warning function (prints to STDOUT in YELLOW/ORANGE).
  echo -e >&1 "${yellow}${underline}WARNING:${rmunderline} ${1:-$(</dev/stdin)}${default}"
}

# Script requires bwa - check for existence:

command -v bwa >/dev/null 2>&1 || { err "BWA doesn't appear to be installed. Aborting."; exit 1; }
command -v samtools >/dev/null 2>&1 || { err "SAMTOOLS doesn't appear to be installed. Aborting."; exit 1; }

# Capture inputs

usage()
{
cat << EOF
usage: $0 options

This script aligns fastq reads to a given
reference, and outputs alignment maps etc.

OPTIONS:
   -h | --help     Show this message
   -r | --ref      Reference fasta sequence
   -1 | --read1    First set of reads ("R1")
   -2 | --read2    Second set of reads ("R2")

EOF
}

# Tolerate long arguments
for arg in "$@" ; do
  shift
  case "$arg" in
    "--help")  set -- "$@" "-h"   ;;
    "--ref")   set -- "$@" "-r"   ;;
    "--read1") set -- "$@" "-1"   ;;
    "--read2") set -- "$@" "-2"   ;;
    *)         set -- "$@" "$arg" ;;
  esac
done
# getopts assigns the arguments to variables
while getopts "hr:1:2:" OPTION ; do
  case $OPTION in
    r) reference=$OPTARG ;;
    1) read1=$OPTARG     ;;
    2) read2=$OPTARG     ;;
    h) usage ; exit      ;;
  esac
done

if [[ -z $reference ]] ; then
  usage
  err "No Reference sequence provided. Exiting." ; exit 1
fi

if [[ -z $read1 ]] ; then
  usage
  err "Forward reads not supplied. Exiting." ; exit 1
fi

if [[ -z $read2 ]] ; then
  usage
  err "Reverse reads not supplied. Exiting." ; exit 1
fi

# Index the reference sequence
bwa index "$reference"
log "BWA referenced index"

log "Running bwa-mem as follows:"
log " -> bwa mem $reference $read1 $read2 > ${reference%.*}_aligned.sam"
bwa mem "$reference"  "$read1" "$read2" > "${reference%.*}"_aligned.sam

log "Created SAM file."

samtools view -S "${reference%.*}"_aligned.sam -b -o ./"${reference%.*}"_aligned.bam
log "Converted SAM to BAM"

samtools sort "${reference%.*}"_aligned.bam "${reference%.*}"_aligned_sorted
log "Sorted BAM file."

samtools index "${reference%.*}"_aligned_sorted.bam
log "Indexed sorted BAM file."


