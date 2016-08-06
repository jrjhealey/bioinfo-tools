#!/bin/bash

# A small script to generate artemis comparison files (nucleic acid comparison)
# since all the webservers are apparently defunct!


# Step 1: Make a BLAST database of the reference sequence:
echo "Please specify the reference sequence (in FASTA format):"
read -e ref_seq
echo "Please provide the name of the BLAST database (also used for the output filenames)":
read -e dbname
echo "Please provide the query sequence to be used in the BLAST step:"
read -e query_seq
echo "Lastly, provide the comparison filename for use with ACT:"
read -e comparison_file

echo -en Executing command with arugments: '\n' Reference - $ref_seq '\n' Database - $dbname '\n' Query sequence - $query_seq '\n' Comparison file - $comparison_file

sleep 3
makeblastdb -in $ref_seq -dbtype 'nucl' -input_type fasta -title $dbname -out $dbname
wait
echo "Database created."

# Step 2: Perform the all-vs-all BLAST using the query sequence and reference database.
wait
blastall -p blastn -d $dbname -i $query_seq -m 8 -q -3 -o $comparison_file
wait
echo "All finished!"

