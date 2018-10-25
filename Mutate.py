#!/usr/bin/env python
from __future__ import print_function
from Bio import SeqIO
from Bio.Seq import Seq
from collections import OrderedDict
import sys
import argparse

# TODO:
#  - create some logic to 'group' mutations that will be applied to the same sequence, to
#    make all switches at once
#    - This will also probably break the verbose transversion output so the maths will need replacing
#  - Create the ability to support INDELS (will also require pairwise alignment so that
#    hamming distances remain meaningful.

def get_args():
    """Parse command line arguments"""

    try:
        parser = argparse.ArgumentParser(
            description='Mutate fasta sequences based on a file of mappings.')
        parser.add_argument('mutation_file', action='store',
                            help='File of mutation mappings like so: "SeqID,X123Y"')
        parser.add_argument('sequences', action='store',
                            help='File of sequences to be mutated (fasta only).')
        parser.add_argument('-v','--verbose', action='store_true',
                            help='Verbose behaviour, printing parameters of the script.')
        parser.add_argument('-o', '--outfile', action='store',
                            help='Output file for mutated sequences (default STDOUT).')
        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            exit(1)
    except:
        sys.stderr.write("An exception occurred with argument parsing. Check your provided options.")

    return parser.parse_args()


class Mutation(object):
    """A class wrapper for sequence IDs so that duplicate IDs can be used in a dictionary"""

    def __init__(self, name):
        self.name = name

    def __repr__(self):
        return "'"+self.name+"'"

    def __str__(self):
        return self.name


def parse_mapfile(mapfile):
    """Return a dict of mapped mutations.
    
    File should resemble:
    
     SequenceID,A123B
     SequenceID2,X234Y
    
    Sequence IDs should exactly match the fasta headers, as parsed by BioPython.
    (">" symbols are optional)
    """

    with open(mapfile, 'r') as handle:
        mut_dict = OrderedDict()
        for line in handle:
            id, change = line.lstrip('>').rstrip('\n').split(',')
            mut_dict[Mutation(id)] = change

    for k,v in mut_dict.items():
        assert v[0].isalpha(), "First character of mutation map is not a valid letter. Got: %s" % v[0]
        assert v[-1].isalpha(), "Last character of mutation map is not a valid letter. Got: %s" % v[-1] 
        assert v[1:-1].isdigit(), "Location string of mutation map is not a valid number. Got: %s" % v[1:-1]

    return mut_dict


def morph(orig, loc, new, mutableseq, verbose):
    """Perform actual sequence change (polymorphism only at present)"""
    # Shift location to offset 0-based index
    loc = loc - 1
    assert mutableseq[loc] == orig, "Sequence does not match the mutation file for pre-exising residue. Expected %s , got %s " % (orig, mutableseq[loc]) 
    if verbose is True:
        print("Performing change: %s -> %s, at location: %d (0 based)" %(orig, new, loc) )
    mutableseq[loc] = new
    return mutableseq


def hamming_distance(s1, s2):
    """Return the Hamming distance between equal-length sequences"""
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))    


def main():
    args = get_args()
    if args.outfile is not None:
        ofh = open(args.outfile, 'w')
    
    # Parse the mutation file (get mutations by sequence)
    mutations = parse_mapfile(args.mutation_file)
    if args.verbose is True:
        print("Got mutations:")
        print(mutations)    
    # Iterate all sequences and make any substitutions necessary
    for record in SeqIO.parse(args.sequences, 'fasta'):
        mutable = record.seq.tomutable()
        for k, v in mutations.items():
            if k.name == record.id:
                orig = v[0]
                new = v[-1]
                loc = int(v[1:-1])
                if args.verbose: print(record.id)               
                newseq = morph(orig, loc, new, mutable, args.verbose)
                if args.verbose is True:
                    print("Original: " + record.seq)
                    print(str((' ' * int(loc-2+11))) + 'V') # Padded string to show where switch happened (not sure how it'll deal with line wrapping 
                    print("New:      " + newseq)
                    print("Distance: " + str(hamming_distance(str(record.seq), str(newseq))))
                if args.outfile is not None:
                    ofh.write(">%s_%s\n%s\n" % (record.id, v, newseq))
                if (args.verbose is False):
                    print(">%s_%s\n%s\n" % (record.id, v, newseq))

    if args.outfile is not None:
        ofh.close()
    
if __name__ == '__main__':
    main()
