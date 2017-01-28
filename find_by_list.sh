#!/bin/bash

# Script that looks for files based on a list in a key file and moves them to a selected location.

# Usage:
# bash find_by_list.sh [listfile] [directory_to_search] [directory_to_copy_to]

# Capture inputs

usage()
{
cat << EOF
usage: $0 options

This script implements "find" to search through a keyfile (a list of 
strings) and copies them to a particular directory.

Implementations:
bash find_by_list.sh -k [keyfile] -p [directory_to_search] -o [directory_to_copy_to]
bash find_by_list.sh --key [keyfile] --path [directory_to_search] --out [directory_to_copy_to]

OPTIONS:
   -h | --help     	  	Show this message
   -k | --key		  	The file of keys/strings to look for.
   						(Each entry of the keyfile must be on a new line)
   -p | --path 	  	  	The top most directory of the tree you
   						wish to search for your queries
   -o | --out			Directory to copy the found directories to
   -f | --file 			Switch find to search for files not directories (default off)	

EOF
}

# Tolerate long arguments
for arg in "$@"; do
  shift
  case "$arg" in
  		"--help")
			set -- "$@" "-h"
			;;
	    	"--key")
    			set -- "$@" "-k"
			;;
   		"--path")
    			set -- "$@" "-p"
    			;;
    		"--file")
			set -- "$@" "-f"
			;;
    		"--out")
       			set -- "$@" "-o"
       			;;
		*)
			set -- "$@" "$arg"
  esac
done
# getopts assigns the arguments to variables
while getopts "hk:p:o:f" OPTION ; do
	case $OPTION in
		k)
		keyfile=$OPTARG
		;;
		p)
		path=$OPTARG
		;;
		o)
		out=$OPTARG
		;;
		f)
		type=$OPTARG
		;;
		h)
		usage
		exit
		;;
	esac
done

if [[ -z $keyfile ]] ; then
    	usage
    	echo "No keyfile was provided. Exiting." ; exit 1
fi
if [[ -z $path ]] ; then
    	usage
    	echo "Path to search was not specified. Exiting." ; exit 1
fi
if [[ -z $out ]] ; then
    	echo "Output directory not supplied. Defaulting to current."
    	out=$(pwd)
    	echo "Directory is ${out}"
fi

# Start the search process
if [[ "$type" == "f"]] ; then
		while read key
			do
        		find "$path" -type f -name "*${key}*" -exec cp -r {} "$out" \;
			done < $keyfile
	else
		if [[ "$type" == "d" ]] || [[ -z "$type" ]]
			then
				while read key
					do
        				find "$path" -type d -name "*${key}*" -exec cp -r {} "$out" \;
				done < $keyfile
		fi
fi


