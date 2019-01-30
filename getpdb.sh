#!/bin/bash
set -eo pipefail

# Script to retrieve PDBs via the command line from the PDB HTTP/FTP

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

usage()
{
cat << EOF
usage: $0 options

This script retrieves PDB files from the Protein DataBank.

OPTIONS:
   -h | --help      Show this message
   -i | --ID        4-letter alphanumeric PDB ID(s) to fetch (e.g. 1abc)
                    To download multiple structures, issue -i each time:
                    ...  -i 1abc -i 2def -i 3ghi etc...
   -m | --mode      Which protocol to use to fetch with ("HTTP" or "FTP")
   -o | --outdir    Directory to move the PDB to (optional, default = ./):
                    $(pwd)
EOF
}

# Tolerate long arguments
for arg in "$@"; do
  shift
  case "$arg" in
   "--help")    set -- "$@" "-h"   ;;
   "--ID")      set -- "$@" "-i"   ;;
   "--mode")    set -- "$@" "-m"   ;;
   "--outdir")  set -- "$@" "-o"   ;;
   *)           set -- "$@" "$arg" ;;
  esac
done
# getopts assigns the arguments to variables
while getopts "hi:m:o:" OPTION ; do
  case $OPTION in
   i) ID+=($OPTARG)   ;;  # appends an array rather than set a single variable
   m) mode=$OPTARG    ;;
   o) outdir=$OPTARG  ;;
   h) usage ; exit 1 ;;
  esac
done
shift $((OPTIND -1))

# Check var and set defaults if necessary


if [[ -z $ID ]] ; then
 usage
 err "No ID was provided. Exiting." ; exit 1
fi

if [[ -z "$mode" ]] ; then
 mode="HTTP"
fi

if [[ -z "$outdir" ]] ; then
 outdir=$(pwd)
 log "Saving in to: ${outdir}"
fi

log "Fetching the following IDs:"
log "${ID[*]}"

for i in "${ID[@]}" ; do
log "Now working on ${i}..."
 if [[ "$mode" == "FTP" ]]
  then
      i=$(echo "$i" |tr '[:upper:]' '[:lower:]')
      wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb${i}.ent.gz 2>&1 | sed '2,$s/^/\      /g' | log ||
        { err "Failed to find a matching PDB entry. Exiting." ; exit 1; }
      gunzip pdb${i}.ent.gz &&  mv -v pdb${i}.ent ${outdir%./}/${i}.pdb | sed '2,$s/^/\         /g' | log
  else
 if [[ "$mode" == "HTTP" ]]
  then
      wget https://files.rcsb.org/download/${i}.pdb 2>&1 | sed '2,$s/^/\      /g' | log ||
        { err "Failed to find a matching PDB entry. Exiting" ; exit 1; }
  fi
 fi
done
