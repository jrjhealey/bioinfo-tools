#!/bin/bash

topdir=$(pwd)

for dir in */ ; do
	cd $dir
	cd nuc/
	less ${dir%/}.clust | grep Aligned | cut -f6 -d' ' | awk -v locus=${dir%/} '{ sum += $1 } ; END { print locus, "Sum of ID:", sum, " Avg ID: ", sum / FNR}' > ${dir%/}_locus_avg_ID.txt
	cd $topdir
done
