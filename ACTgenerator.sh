#!/bin/bash
# A small script to generate artemis comparison files (nucleic acid comparison)
# since all the webservers are apparently defunct!
set -eo pipefail

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
  # Logging function.
  # Prints to STDOUT in WHITE
  echo -e >&1 "${white}${underline}INFO:${rmunderline} $1${default}"
}

err(){
  # Error function
  # Prints to STDERR in RED
  echo -e >&2 "${red}${underline}ERROR:${rmunderline} $1${default}"
}

warn(){
  # Warning function
  # Prints to STDOUT in YELLOW/ORANGE
  echo -e >&1 "${yellow}${underline}WARNING:${rmunderline} $1${default}"
}

usage(){
  # Capture inputs
cat << EOF >&2
usage: $0 options

This script generates the necessary BLAST comparison file
for use with the Artemis Comparison tool, when comparing 2 genomes.

A typical invocation might look like:

 $ bash ACTgenerator.sh -r reference.fasta -q query.fasta -d databasename -o outputdir -t

OPTIONS:
   -h | --help      Show this message.
   -r | --ref       Reference fasta sequence.
   -q | --query     Query fasta sequence.
   -d | --database  Name of the BLAST database comparison file.
   -o | --outdir    Directory to output all the files to.
   -t | --tidy      Make the script tidy up after itself (T/F).
                    ${red}(Default True.)${default}

This script will take 2 input sequences (in FASTA at present) and creates the BLAST
comparison file that is used by the ARTEMIS Comparison tool (ACT). At the moment, it
only performs a 1-vs-1 comparison, but may be expanded in the future.

By default, the script requires only input query and reference sequences, with all
other options taking on default values. The script will also remove any intermediate
files by default (that includes if not specified. To retain these files provide the -t
flag with a "False" argument (case insensitive).
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

if [[ $tidy =~ ^[Tt][Rr][Uu][Ee]$ ]] || [[ $tidy =~ ^[Tt]$ ]] || [[ -z $tidy ]] ; then
  warn "Tidying database files from ${outdir}."
  rm -v "${outdir%/}"/"${database}"*
elif [[ $tidy =~ ^[Ff][Aa][Ll][Ss][Ee]$ ]] || [[ $tidy =~ ^[Ff]$ ]] ; then
  :   # Do nothing if false-y
else
 err 'Unrecognised argument to tidy (should be T(rue), F(alse) or empty). Leaving files unmodified.'
fi

