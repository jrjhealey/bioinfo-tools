#!/bin/bash
set -eo pipefail
tr=$(tput setaf 1) # text = red
ty=$(tput setaf 3)
tw=$(tput setaf 7)
smul=$(tput smul)
rmul=$(tput rmul)
df=$(tput sgr0)    # text = reset
# A small script to generate artemis comparison files (nucleic acid comparison)
# since all the webservers are apparently defunct!

log(){
  # Logging function.
  # Prints to STDOUT in WHITE
  echo -e >&1 "${tw}${smul}INFO:${rmul} $1${df}"
}

err(){
  # Error function
  # Prints to STDERR in RED
  echo -e >&2 "${tr}${smul}ERROR:${rmul} $1${df}"
}

warn(){
  # Warning function
  # Prints to STDOUT in YELLOW/ORANGE
  echo -e >&1 "${ty}${smul}WARNING:${rmul} $1${df}"
}

usage(){
  # Capture inputs
cat << EOF >&2
usage: $0 options

This script generates the necessary BLAST comparison file
for use with the Artemis Comparison tool, when comparing 2 genomes.

A typical invocation might look like:

 $ bash ACTgenerator.sh -r reference.fasta -q query.fasta -d databasename -o outputdir -t "T"

OPTIONS:
   -h | --help      Show this message.
   -r | --ref       Reference fasta sequence.
   -q | --query     Query fasta sequence.
   -d | --database  Name of the BLAST database comparison file.
   -o | --outdir    Directory to output all the files to.
   -t | --tidy      Make the script tidy up after itself (T/F).
                    ${tr}Default F.${df}

EOF
}

# Tolerate long arguments
for arg in "$@"; do
 shift
 case "$arg" in
   "--help")      set -- "$@" "-h"   ;;
   "--reference") set -- "$@" "-r"   ;;
   "--query")     set -- "$@" "-q"   ;;
   "--database")  set -- "$@" "-d"   ;;
   "--outdir")    set -- "$@" "-o"   ;;
   "--tidy")      set -- "$@" "-t"   ;;
   *)             set -- "$@" "$arg" ;;
 esac
done
# getopts assigns the arguments to variables
while getopts "hr:q:d:o:t:" OPTION ;do
  case $OPTION in
    r) reference=$OPTARG ;;
    q) query=$OPTARG     ;;
    d) database=$OPTARG  ;;
    o) outdir=$OPTARG    ;;
    t) tidy=$OPTARG      ;;
    h) usage; exit 0     ;;
	esac
done

# If no args, show help
if [[ $# -eq 0 ]] ; then
    usage ; exit 1
fi

if [[ -z $reference ]]; then
 usage ; err "No Reference sequence provided. Exiting." ; exit 1
fi

if [[ -z $query ]]; then
 usage ;  err "Query not supplied. Exiting." ; exit 1
fi

if [[ -z $database ]]; then
 database="${reference%.*}_DB"
 warn "No database name was specified. Using the reference sequence name and appending _DB."
fi

if [[ -z $outdir ]]; then
 outdir=$(pwd)
 warn "No output directory was specified. Defaulting to the current working directory: $outdir"
fi

#####

# Step 1: Make a BLAST database of the reference sequence:
# Script requires blastn (NOT LEGACY BLAST) and makeblastdb in path - check for existence:
command -v makeblastdb >/dev/null 2>&1 || { err "makeblastdb doesn't appear to be installed. Aborting."; exit 1; }
command -v blastn >/dev/null 2>&1 || { err "BLAST+ doesn't appear to be installed. Aborting."; exit 1; }


log "Running makeblastdb:"
log " -> makeblastdb -in "$reference" -dbtype 'nucl' -title "$database" -out "${outdir%/}"/"${database}" -parse_seqids"

makeblastdb -in "$reference" -dbtype 'nucl' -title "$database" -out "${outdir%/}"/"${database}" -parse_seqids
log "Database created."

# Step 2: Perform the all-vs-all BLAST using the query sequence and reference database.
log " -> blastn -db "${outdir%/}"/"${database}" -query "$query" -outfmt 6 -out "${outdir%/}"/"${query%.*}"_vs_"${reference%.*}".act"
blastn -db "${outdir%/}"/"${database}" -query "$query" -outfmt 6 -out "${outdir%/}"/"${query%.*}"_vs_"${reference%.*}".act

log "All finished! The comparison file is called: ${query%.*}_vs_${reference%/*}.act"

if [[ $tidy =~ ^[Tt][Rr][Uu][Ee]$ ]] || [[ $tidy =~ ^[Tt]$ ]] ; then
  warn "Tidying database files from ${outdir}."
  rm -v "${outdir%/}"/"${database}"*
elif [[ $tidy =~ ^[Ff][Aa][Ll][Ss][Ee]$ ]] || [[ $tidy =~ ^[Ff]$ ]] || [[ -z $tidy ]] ; then
  :   # Do nothing if false-y
else
 err 'Unrecognised argument to tidy (should be T(rue), F(alse) or empty). Leaving files unmodified.'
fi

