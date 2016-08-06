#!/bin/bash

# Script to prepare folders and files for itasser
dir="$1"

for file in $(ls "$dir") ; do
        mkdir "${dir%/}"/"${file%.*}"
        echo File is $file
        mv -v "${dir%/}"/"$file" "${dir%/}"/"${file%.*}"/seq.fasta
done
