#!/bin/bash

paste $1 $2 |\
awk '{ printf("%s",$0); n++; if(n%4==0) { printf("\n");} else { printf("\t\t");} }' |\
shuf  |\
head |\
sed 's/\t\t/\n/g' |\
awk -F '\t' '{print $1 > "file1.fastq"; print $2 > "file2.fatsq"}'
