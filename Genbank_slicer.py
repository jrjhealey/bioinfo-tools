#!/usr/bin/env python

# This script is designed to take a genbank file and 'slice out'/'subset'
# regions (genes/operons etc.) and produce a separate file. This can be
# done explicitly by telling the script which base sites to use, or can
# 'decide' for itself by blasting a fasta of the sequence you're inter-
# ed in against the Genbank you want to slice a record out of.

# Note, the script (obviously) does not preseve the index number of the
# bases from the original

# Based upon the tutorial at:
# http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc44

# This script depends on BLAST and having the makeblastdb functionality
# available if BLAST_MODE is active. It also depends on Biopython.

# Set up and handle arguments:

import warnings
import sys
import subprocess
import os
import time

try:
    from Bio import SeqIO
except ImportError:
    msg = """
Could not import the BioPython module which means it is probably
not installed, or at least not available in the PYTHONPATH for this 
particular binary.

If you have conda (recommended) try running:

 $ conda install -c anaconda biopython
 
or, alternatively with python/pip:

 $ python -m pip install biopython
"""
    sys.stderr.write(msg)
    sys.exit(1)

__author__ = "Joe R. J. Healey"
__version__ = "1.2.0"
__title__ = "Genbank_slicer"
__license__ = "GPLv3"
__author_email__ = "J.R.J.Healey@warwick.ac.uk"

# TODO:
# - Alter the script to slice other sequence types by removing the
#   requirement for Genbanks.
# - Figure out how to deal with multigenbanks

# Import SearchIO and suppress experimental warning
from Bio import BiopythonExperimentalWarning

with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonExperimentalWarning)
    from Bio import SearchIO


def convert(basename, genbank):
    """Convert the provided genbank to a fasta to BLAST."""

    refFasta = "{}.fasta.tmp".format(basename)
    SeqIO.convert(genbank, "genbank", refFasta, "fasta")

    return refFasta


def runBlast(basename, refFasta, fasta, verbose):
    """Synthesise BLAST commands and make system calls"""

    resultHandle = "{}.blastout.tmp".format(basename)
    blastdb_cmd = "makeblastdb -in {0} -dbtype nucl -title temp_blastdb".format(
        refFasta
    )
    blastn_cmd = "blastn -query {0} -strand both -task blastn -db {1} -perc_identity 100 -outfmt 6 -out {2} -max_target_seqs 1".format(
        fasta, refFasta, resultHandle
    )

    print("Constructing BLAST Database: " + "\n" + blastdb_cmd)
    print("BLASTing: " + "\n" + blastn_cmd)
    DB_process = subprocess.Popen(
        blastdb_cmd,
        shell=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    DB_process.wait()
    BLAST_process = subprocess.Popen(
        blastn_cmd,
        shell=True,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    BLAST_process.wait()

    return resultHandle


def getIndices(resultHandle):
    """If not provided directly by the user, this function retrieves the best BLAST hit's indices."""

    blast_result = SearchIO.read(resultHandle, "blast-tab")

    print(blast_result[0][0])
    start = blast_result[0][0].hit_start
    end = blast_result[0][0].hit_end

    return start, end


def slice(start, end, genbank, FPoffset, TPoffset):
    """Subset the provided genbank to return the sub record."""

    try:
        seqObj = SeqIO.read(genbank, "genbank")
    except ValueError:
        sys.stderr.write(
            "There is more than one sequence in the target sequence file.\n"
            "This script requires that there be only 1 currently, else the retrieved indices are meaningless.\n"
            "Please concatenate the target sequence and try again."
            "The script will carry on with just the first sequence record."
        )
        seqObj = list(SeqIO.parse(genbank, "genbank"))[0]

    subRecord = seqObj[start - FPoffset : end + TPoffset]

    return subRecord


def main():
    ###################################################################################################
    # Parse command line arguments
    import argparse

    try:
        parser = argparse.ArgumentParser(
            description="This script slices entries such as genes or operons out of a genbank, subsetting them as their own file."
        )
        parser.add_argument(
            "-g",
            "--genbank",
            action="store",
            required=True,
            help="The genbank file you wish to subset.",
        )
        parser.add_argument(
            "-o",
            "--outfile",
            action="store",
            help="If specifed, the script will write a file, otherwise redirect STDOUT for pipelines.",
        )
        parser.add_argument(
            "-s", "--start", type=int, help="Integer base index to slice from."
        )
        parser.add_argument(
            "-e", "--end", type=int, default=0, help="Integer base index to slice to."
        )
        parser.add_argument(
            "-F",
            "--FPoffset",
            action="store",
            type=int,
            default=0,
            help="If you want to capture regions around the target site, specify the 5' offset.",
        )
        parser.add_argument(
            "-T",
            "--TPoffset",
            action="store",
            type=int,
            default=0,
            help="If you want to capture regions around the target site, specify the 3' offset.",
        )
        parser.add_argument(
            "-b",
            "--blastmode",
            action="store_true",
            help="If flag is set, provide an input fasta (-f | --fasta), and BLAST will be called to find your sequence indices in the original genbank.",
        )
        parser.add_argument(
            "-f",
            "--fasta",
            action="store",
            help="The operon fasta to pull annotations from the provided genbank.",
        )
        parser.add_argument(
            "-t",
            "--tidy",
            action="store_true",
            help='Tell the script whether or not to remove the temporary files it generated during processing. Off by default. WARNING: removes files based on the "tmp" string.',
        )
        parser.add_argument(
            "-m",
            "--meta",
            action="store",
            default=None,
            help="Specify a string to use as the Genbank meta-data fields if the parent one doesn't contain anything. Else inherits from parent sequence.",
        )
        parser.add_argument(
            "-v",
            "--verbose",
            action="store_true",
            help="Verbose behaviour. Shows progress of BLAST etc. if appropriate.",
        )

        args = parser.parse_args()

    except NameError:
        sys.stderr.write(
            "An exception occured with argument parsing. Check your provided options."
        )

    genbank = args.genbank
    fasta = args.fasta
    split = os.path.splitext(args.genbank)
    basename = os.path.basename(split[0])
    start = args.start
    end = args.end
    FPoffset = args.FPoffset
    TPoffset = args.TPoffset
    blastMode = args.blastmode
    outfile = args.outfile
    tidy = args.tidy
    meta = args.meta
    verbose = args.verbose

    # Main code:
    if (blastMode is not False) and (fasta is not None):
        refFasta = convert(basename, genbank)
        resultHandle = runBlast(basename, refFasta, fasta, verbose)
        start, end = getIndices(resultHandle)
    elif (blastMode is not False) and (fasta is None):
        sys.stderr.write("No fasta was provided so BLAST mode cannot be used.")
        sys.exit(1)

    if start is None or end is None:
        sys.stderr.write(
            "No slice indices have been specified or retrieved from blastout."
        )
        sys.exit(1)

    subRecord = slice(start, end, genbank, FPoffset, TPoffset)
    # Populate the metadata of the genbank
    if meta is not None:
        subRecord.id = meta
        subRecord.locus = meta + "subregion"
        subRecord.accession = meta
        subRecord.name = meta
        subRecord.annotations["date"] = time.strftime("%d-%b-%Y").upper()
    # Other field options include:
    #  subRecord.annotations["source"]
    # 			["taxonomy"]
    # 			["keywords"]
    # 			["references"]
    # 			["accessions"]
    # 			["data_file_division"] e.g. BCT, PLN, UNK
    # 			["organism"]
    # 			["topology"]

    if outfile is not None:
        directory = os.path.dirname(outfile)
        if not os.path.exists(directory):
            os.makedirs(directory)

        SeqIO.write(subRecord, outfile, "genbank")
    else:
        print(subRecord.format("genbank"))

    if blastMode is True and tidy is True:
        if verbose is True:
            rm = "rm -v ./*tmp*"
        else:
            rm = "rm ./*tmp*"

        subprocess.Popen(rm, shell=True)


if __name__ == "__main__":
    main()
