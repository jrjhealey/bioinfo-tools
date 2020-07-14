#!/usr/bin/env python

"""
Simple script that uses Biopython AlignIO to convert
between various file formats for sequence alignments.

See documentation online at:
http://biopython.org/wiki/AlignIO
"""


# This program is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as
# published by the Free Software Foundation, either version 3 of
# the License, or (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.


import traceback
import os
import sys
import warnings

__author__ = "Joe R. J. Healey"
__version__ = "1.0.1"
__title__ = "x2phy"
__license__ = "GPLv3"
__author_email__ = "jrj.healey@gmail.com"


def parseArgs():
    """Parse command line arguments"""

    import argparse
    try:
        parser = argparse.ArgumentParser(description='Convert between alignment file types supported by BioPython')
        parser.add_argument('-i',
    	                    '--infile',
                      	    action='store',
                            default='None',
                            help='The input phylogeny to convert.')
        parser.add_argument('-j',
                            '--intype',
                            action='store',
			    default='None',
                            help='Specify the input file type. Script will try to guess from the extension if it is standard.')
        parser.add_argument('-o',
                            '--outfile',
                            action='store',
		            default='None',
                            help='The output file name (default is infile with new extension.)')
        parser.add_argument('-p',
                            '--outtype',
                            action='store',
                            help='The output file type.')
        parser.add_argument('-v',
                           '--verbose',
                           action='count',
                           default=0,
                           help='Verbose behaviour, emit messages about progress/debugging etc.')
        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            exit(1)

    except NameError:
        sys.stderr.write("An exception occurred with argument parsing. Check your provided options.")
        traceback.print_exc()

    return parser.parse_args()


def convert(infile, type, outtype, outfile):
    """Make AlignIO call to convert using the specified parameters"""

    from Bio import AlignIO

    ifh = AlignIO.parse(infile, type)
    AlignIO.write(ifh, outfile, outtype)


def getOutfile(infile, outtype):
    """Synthesise a default output name if none was provided"""
    return os.path.splitext(infile)[0] + "." + outtype


def guessExt(infile, verbose):
    """If no input type was specified, guess it from the extension name or emit a warning"""
    extension = os.path.splitext(infile)[1]
    if verbose > 0:
        sys.stderr.write("Extension is " + extension)
    # Figure out what extension to return
    if extension in (".clust", ".clustal",".aln"):
        type = "clustal"
    elif extension in (".emb", ".emboss"):
        type = "emboss"
    elif extension in (".fa", ".fasta", ".fas", ".fna", ".faa", ".afasta"):
        type = "fasta"
    elif extension in (".phy", ".phylip"):
        type = "phylip"
    elif extension in (".nexus", ".paup"):
        type = "nexus"
    else:
        sys.stderr.write("Couldn't determine the file type from the extension. Reattempt with the -j|--intype option specified.")
        sys.exit(1)

    if verbose > 0:
        sys.stderr.write("Your file looks like a " + type)
    return type



def main():
    args = parseArgs()
    if args.verbose > 1:
        print(args)

    if args.intype is 'None':
        args.intype = guessExt(args.infile, args.verbose)
    if args.outfile is 'None':
        args.outfile = getOutfile(args.infile, args.outtype)

    if args.verbose > 0:
        sys.stderr.write("Converting " + args.infile + " of type: " + args.intype + ", to " + args.outfile + " of type: " + args.outtype)
    if args.infile is not 'None':
        convert(args.infile, args.intype, args.outtype, args.outfile)
    else:
        sys.stderr.write("No input file was provided. Check the options you provided.")
        sys.exit(1)

if __name__ == '__main__':
    main()
