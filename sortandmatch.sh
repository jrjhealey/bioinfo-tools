#!/bin/bash

# A small script to take 2 files and sort File 2 by the order of strings in File 1.
# It's behaviour is somewhat similar to that of VLOOKUP.

usage()
{
cat << EOF
usage: $0 options

A small script to take 2 files and sort File 2 by the order of strings in File 1, and
then join the 2 files together on a row-row basis.
A particular column in File 1 to sort by can be specified.

Behaviour is somewhat similar to that of Excel's VLOOKUP function.

OPTIONS:
   -h | --help     	Show this message
   -k | --keyfile	A file specifying the desired order (a column of strings
			which will be common to both files).
   -i | --infile	The file you wish to sort.
   -f | --field		The field you wish to sort based upon in the keyfile [def = 1].
EOF
}

# Tolerate long arguments
for arg in "$@"; do
  shift
  case "$arg" in
  		"--help")
		set -- "$@" "-h"
		;;
	    	"--keyfile")
    		set -- "$@" "-k"
		;;
    		"--infile")
       		set -- "$@" "-i"
       		;;
		"--field")
		set -- "$@" "-f"
		;;
	    	*)
		set -- "$@" "$arg"
  esac
done
# getopts assigns the arguments to variables
while getopts "hk:i:f:" OPTION
do
	case $OPTION in
		k)
		keyfile=$OPTARG
		;;
		i)
		infile=$OPTARG
		;;
		f)
		indexfield=$OPTARG
		;;
		h)
		usage
		exit 1
		;;
	esac
done

if [[ -z $keyfile ]] ; then
    	usage
    	echo "No keyfile provided. Exiting." ; exit 1
fi
if [[ -z $infile ]] ; then
    	usage
    	echo "Infile not supplied. Exiting." ; exit 1
fi
if [[ -z $indexfield ]] ; then
    	indexfield=1
fi

# sort the file with awk, directed in to a paste command which unites the sorted
# File2 and the keyfile line-by-line.

# Cut --complement retains all columns except the specified one, removing the duplication
# effect of pasteing the 2 files that share a column together.

paste "$keyfile" <(awk -v field="$indexfield" 'NR==FNR{o[FNR]=$field; next} {t[$1]=$0} END{for(x=1; x<=FNR; x++){y=o[x]; print t[y]}}' "$keyfile" "$infile") | cut --complement -f "$indexfield"


