#!/usr/bin/python

import os
import sys

inFile = open(sys.argv[1],'r')

totalBP = 0
gcBP = 0
headerCount = 0


for line in inFile:
    if line[0] == ">":
        headerCount += 1
    else:
        seqLine = line.strip().lower()
        totalBP += len(seqLine)
        gcBP += seqLine.count('g')
        gcBP += seqLine.count('c')

print ("Locus %s" % str(sys.argv[1])) + ' GC: ' + str(float(gcBP) / totalBP)
