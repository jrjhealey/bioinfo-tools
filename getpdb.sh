#!/bin/bash

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

# Tolerate long arguments
for arg in "$@"; do
  shift
  case "$arg" in
  		"--help")
		set -- "$@" "-h"
		;;
	    	"--ID")
    		set -- "$@" "-i"
		;;
   		"--mode")
    		set -- "$@" "-m"
    		;;
		"--outdir")
		set -- "$@" "-o"
		;;
	    	*)
		set -- "$@" "$arg"
  esac
done
# getopts assigns the arguments to variables
while getopts "hi:m:o:" OPTION ; do
	case $OPTION in
		i)
		ID+=($OPTARG)
		;;
		m)
		mode=$OPTARG
		;;
		o)
		outdir=$OPTARG
		;;
		h)
		usage
		exit
		;;
	esac
done
shift $((OPTIND -1))

if [[ -z $ID ]]
	then
    	usage
    	echo "No ID was provided. Exiting." ; exit 1
fi
if [[ -z "$mode" ]]
	then
 	mode="HTTP"
fi
if [[ -z "$outdir" ]]
	then
	outdir=$(pwd)
	echo "Saving in to ${outdir}"
fi

echo "The first value of array 'IDs' is $ID"
echo "The last value of array 'IDs' is ${ID[-1]}"
echo "The entirety of array 'IDs' is ${ID[@]}"
#####

# FTP fetch
for i in "${ID[@]}" ; do
 if [[ "$mode" == "FTP" ]]
  then
      wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb${i}.ent.gz
      gunzip pdb${i}.ent.gz &&  mv -v pdb${i}.ent ${outdir%./}/${i}.pdb
  else
 if [[ "$mode" == "HTTP" ]]
  then
      wget http://www.rcsb.org/pdb/files/${i}.pdb.gz
      gunzip -c ${i}.pdb.gz > "${outdir%./}"/${i}.pdb && rm ${i}.pdb.gz
  fi
 fi
done
