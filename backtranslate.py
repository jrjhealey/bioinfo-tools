#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Py3 and UTF-8 encoding required for special characters.

# Code that converts a protein sequence in to a representation of the available
# codon encodings. Note, it doesnt enumerate all individual possibilities (as one
# might do using itertools.product as this simply takes too long to calculate/print).

# Usage:
#  $ python3 backtranslate.py <outputformat> <proteinfile> <proteinfileformat>
# E.g:
#  $ python3 backtranslate.py pretty myproteins.fasta fasta
# If the output is wrapping weirdly:
#  $ python3 backtranslate.py pretty myproteins.fasta fasta |cut -c1-$(stty size </dev/tty | cut -d' ' -f2)

#TODO:
# - Refactor to classes/structures
# - Add argparse.
# - Include different codon tables (create a codon table class to hold them all?)
#   - Maybe lean on Bio.IUPAC for this?
# - Figure out screen wrapping. Currenty needs exporting to a file, or trimming with
#   the construct ( |cut -c1-$(stty size </dev/tty | cut -d' ' -f2) )


# Big thanks to the guys on Code Golf and Programming Puzzle SE for the pretty/boxed output

codon_table_11 = {'A': ['GCU', 'GCC', 'GCA', 'GCG'],
                  'C': ['UGU', 'UGC'],
                  'D': ['GAU', 'GAC'],
                  'E': ['GAA', 'GAG'],
                  'F': ['UUU', 'UUC'],
                  'G': ['GGU', 'GGC', 'GGA', 'GGG'],
                  'H': ['CAU', 'CAC'],
                  'I': ['AUU', 'AUC', 'AUA'],
                  'K': ['AAA', 'AAG'],
                  'L': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
                  'M': ['AUG'],
                  'N': ['AAU', 'AAC'],
                  'P': ['CCU', 'CCC', 'CCA', 'CCG'],
                  'Q': ['CAA', 'CAG'],
                  'R': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
                  'S': ['AGU', 'AGC', 'UCU', 'UCC', 'UCA', 'UCG'],
                  'T': ['ACU', 'ACC', 'ACA', 'ACG'],
                  'V': ['GUU', 'GUC', 'GUA', 'GUG'],
                  'W': ['UGG'],
                  'X': ['nnn'],   # Tolerates unknown AAs but returns unknown codon
                  'Y': ['UAU', 'UAC'],
                  '*': ['UAA', 'UAG', 'UGA']}


def parse(infile, format):
    for protein in SeqIO.parse(infile, format):
        for i in protein:
            amino_acids = []
            for residue in protein:
                amino_acids.append(codon_table_11[residue])
    yield amino_acids, protein.id

def comboprint(matrix):
    m = len(max(matrix, key = len))
    s = [""] * m
    for vector in matrix:
       for i in range(m):
           s[i] += (vector[i] if i < len(vector) else " " * max( [len(x) for x in vector] ) ) + " "
    print('\n'.join(s))

def comboprint_boxed(matrix):
    length = max(map(len, matrix))
    lengths = []
    for index in range(len(matrix)):
        matrix[index] = ([""] * -((len(matrix[index]) - length) // 2) + matrix[index] + ([""] * ((length - len(matrix[index])) // 2)))
        hlength = max(map(len, matrix[index]))
        lengths.append(hlength)
        matrix[index] = [item.ljust(hlength) for item in matrix[index]]
    horiz = ["─" * e for e in lengths]
    print("┬".join(horiz).join("┌┐"))
    for row in list(zip(*matrix))[:-1]:
        print("│".join(row).join("││"))
        print("┼".join(horiz).join("├┤"))
    print("│".join(col[-1] for col in matrix).join("││"))
    print("┴".join(horiz).join("└┘"))

def comboprint_pretty(matrix):
    length = max(map(len, matrix))
    for index in range(len(matrix)):
        matrix[index] = ([""] * -((len(matrix[index]) - length) // 2) + matrix[index] + ([""] * ((length - len(matrix[index])) // 2)))
        hlength = max(map(len, matrix[index]))
        matrix[index] = [item.ljust(hlength) for item in matrix[index]]
    for row in list(zip(*matrix)):
        print(" ".join(row))


if __name__ == "__main__":
    import sys
    from Bio import SeqIO
    if not sys.argv[3]:
        sys.argv[3] = "fasta"

    if sys.version_info < (3, 0):
        sys.exit("Exited (1). Requires python3 to run correctly.")

    functions = {'pretty': comboprint_pretty,
                'boxed': comboprint_boxed,
                'simple': comboprint}

    for amino_acids, name in parse(sys.argv[2], sys.argv[3]):
        print(name)
        functions.get(sys.argv[1], lambda: "Invalid function choice {pretty|boxed|simple}")(amino_acids)
