#!/bin/bash
set -eo pipefail

# Script to retrieve PDBs via the command line from the PDB HTTP/FTP

# Capture inputs

usage()
{
cat << EOF
usage: $0 options

This script retrieves PDB files from the Protein DataBank.

OPTIONS:
   -h | --help             Show this message
   -i | --ID		   4-letter alphanumeric PDB ID(s) to fetch (e.g. 1abc)
			   To download multiple structures, issue -i each time:
			     ...  -i 1abc -i 2def -i 3ghi etc...
   -m | --mode		   Which protocol to use to fetch with ("HTTP" or "FTP")
   -o | --outdir	   Directory to move the PDB to (optional, default = ./):
   			   $(pwd)

EOF
}

log(){
echo -e >&1 "\033[32mINFO: $1\033[0m"
}

err(){
echo -e >&2 "\033[31mERROR: $1\033[0m"
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
   h) usage && exit 1 ;;
  esac
done
shift $((OPTIND -1))

# Check var and set defaults if necessary
if [[ -z $ID ]] ; then
 usage
 echo "No ID was provided. Exiting." ; exit 1
fi

if [[ -z "$mode" ]] ; then
 mode="HTTP"
fi

if [[ -z "$outdir" ]] ; then
 outdir=$(pwd)
 echo "Saving in to ${outdir}"
fi

log "Fetching the following IDs:"
printf '%s\n' "${ID[@]}"
# TODO:
#  Add logic for retrieval from mmCIF?
for i in "${ID[@]}" ; do
log "Now working on ${i}..."
 if [[ "$mode" == "FTP" ]]
  then
      i=$(echo "$i" |tr '[:upper:]' '[:lower:]')
      wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb${i}.ent.gz || { err "Failed to find a matching PDB entry. Exiting." ; exit 1; }
      gunzip pdb${i}.ent.gz &&  mv -v pdb${i}.ent ${outdir%./}/${i}.pdb
  else
 if [[ "$mode" == "HTTP" ]]
  then
      wget https://files.rcsb.org/download/${i}.pdb || { err "Failed to find a matching PDB entry. Exiting" ; exit 1; }
  fi
 fi
done
