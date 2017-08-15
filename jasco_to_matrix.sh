#!/bin/bash

# Takes a JASCO 1.5 format plaintext file, extracts the measurements
# and creates a matrix for downstream plotting and analysis

# ENSURE FILENAMES HAVE NO SPACES IN BEFORE PROCEEDING.

# Usage: bash jasco_to_matrix.sh outfilename spectra1 spectra2 ... spectraN

outfile=$1

jasco=( "${@:2}" )

for file in "${jasco[@]}" ; do
 cat $file | grep -E "^([0-9]+\.[0-9]+)" | cut -f 1-2 > "${file%.txt}"_measurements.txt
done


echo ",20,25,30,35,40,45,50,55,60,65,70,75,80,85,90" > "$outfile"
paste -d '\t' *_measurements.txt | cut -f 1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30 | gsed 's/\t/,/gi' >> "$outfile"
