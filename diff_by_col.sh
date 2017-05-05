#!/bin/bash

# Show lines between 2 files that differ in a given column

file1=$1
file2=$2
indexfield=$3

echo "Lines in file 2 that vary in column $indexfield file 1:"
echo "-------------------------------------------------------"
awk -v field="$indexfield" 'NR==FNR{c[$field]++;next};c[$field] == 0' $file2 $file1
echo "Lines in file 1 that vary in column $indexfield file 2:"
echo "-------------------------------------------------------"
awk -v field="$indexfield" 'NR==FNR{c[$field]++;next};c[$field] == 0' $file1 $file2
