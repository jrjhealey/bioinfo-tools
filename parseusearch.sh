#!/bin/bash

#A script to parse the output blast file from usearch.

# Usage: parseusearch.sh inputblastfile [true] > sortedblastfile.txt
# (If redirect isn't provided, the input file will be sorted).
# if [true] is provided (without brackets), the script will also output
# a file which is just a list of matched locustags.


# Sort results by blast identity to the query

sort -k3 -nr $1 > ${1%.*}_sortedresult.txt

# A list of matched locus tag stripped from the above file.
#if [ -z "$2" ] && [ "$2" == "true" ];
#then
#	cut -f1 $1 > Locus_tags.txt
#else 
#	echo "Variable 2 was either not set or set to false."
#fi
