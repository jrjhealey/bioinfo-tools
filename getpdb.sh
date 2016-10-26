#!/bin/bash

# Script to retrieve PDBs via the command line from the PDB HTTP/FTP

# Capture inputs

usage()
{
cat << EOF
usage: $0 options

This script retrieves PDB files from the Protein DataBank.

OPTIONS:
   -h | --help     	   Show this message
   -i | --ID		   4-letter alphanumeric PDB ID to fetch (e.g. 3izo)
   -m | --mode		   Which protocol to use to fetch with ("HTTP" or "FTP")
   -o | --outdir	   Directory to move the PDB to (optional, default = ./)
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
while getopts "hi:m:o:" OPTION
do
	case $OPTION in
		i)
		ID=$OPTARG
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
#####

# FTP fetch
if [[ "$mode" == "FTP" ]]
 then
      wget ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb${ID}.ent.gz
      gunzip pdb${ID}.ent.gz &&  mv -v pdb${ID}.ent ${outdir%./}/${ID}.pdb
 else
 if [[ "$mode" == "HTTP" ]]
 then 
      wget http://www.rcsb.org/pdb/files/${ID}.pdb.gz
      gunzip -c ${ID}.pdb.gz > "${outdir%./}"/${ID}.pdb && rm ${ID}.pdb.gz
 fi
fi

