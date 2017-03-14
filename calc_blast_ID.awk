#!/bin/awk

{s += $3}
END {print "File analysed was: - " FILENAME ". Sum of ID:", s, " Avg ID: ", s / FNR}

