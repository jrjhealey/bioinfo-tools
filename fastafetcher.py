#!/usr/bin/env python

# Extract fasta files by their descriptors stored in a separate file.
# Requires biopython

# TODO:
# - Create more sophisticated logic for matching IDs/Descriptions/Partial matches etc.
#    - Create a mode variable to encapsulate invert/partial/description/id etc?
from Bio import SeqIO
import sys, six
import argparse


def get_keys(args):
    """Turns the input key file into a list. May be memory intensive."""

    with open(args.keyfile, "r") as kfh:
        keys = [line.rstrip("\n").lstrip(">") for line in kfh]
    return keys


def get_args():
    try:
        parser = argparse.ArgumentParser(
            description="Retrieve one or more fastas from a given multifasta."
        )
        parser.add_argument(
            "-f",
            "--fasta",
            action="store",
            required=True,
            help="The multifasta to search.",
        )
        parser.add_argument(
            "-k",
            "--keyfile",
            action="store",
            help="A file of header strings to search the multifasta for. Must be one per line.",
        )
        parser.add_argument(
            "-s",
            "--string",
            action="store",
            help="Provide a string to look for directly, instead of a file (can accept a comma separated list of strings).",
        )
        parser.add_argument(
            "-o",
            "--outfile",
            action="store",
            default=None,
            help="Output file to store the new fasta sequences in. Just prints to screen by default.",
        )
        parser.add_argument(
            "-v",
            "--verbose",
            action="store_true",
            help="Set whether to print the key list out before the fasta sequences. Useful for debugging.",
        )
        parser.add_argument(
            "-i",
            "--invert",
            action="store_true",
            help="Invert the search, and retrieve all sequences NOT specified in the keyfile.",
        )
        parser.add_argument(
           "-m",
           "--method",
           action="store",
           choices=["exact", "partial"],
           default="exact",
           help="Search the headers as exact matches, or as partial substring matches. "
                "The latter is dangerous, as headers may be matched twice, so be sure "
                "your headers/keys are unique to their respective sequences."
       )

        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            sys.exit(1)

    except NameError:
        sys.stderr.write(
            "An exception occured with argument parsing. Check your provided options."
        )
        sys.exit(1)

    return parser.parse_args()


def main():
    """Takes a string or list of strings in a text file (one per line) and retreives them and their sequences from a provided multifasta."""
    args = get_args()
    # Call getKeys() to create the list of keys from the provided file:
    if not (args.keyfile or args.string):
        sys.stderr.write("No key source provided. Exiting.")
        sys.exit(1)
    if args.keyfile:
        keys = get_keys(args)
    else:
        keys = args.string.split(",")

    if args.verbose:
        if args.invert is False:
            sys.stderr.write("Fetching the following keys:\n")
            for key in keys:
                sys.stderr.write(key + "\n")
        elif args.invert is True:
            sys.stderr.write(
                "Ignoring the following keys, and retreiving everything else from: {}\n".format(
                    args.fasta
                )
            )
            for key in keys:
                sys.stderr.write(key + "\n")
        sys.stderr.write(
            "------------------------------------------------------------\n"
        )

    # Parse in the multifasta and assign an iterable variable:
    to_write = []
    for seq in SeqIO.parse(args.fasta, "fasta"):
        if args.invert is False:
            if args.method == "exact":
                if seq.id in keys:
                    print(seq.format("fasta"))
                    to_write.append(seq)

            elif args.method == "partial":
                if any(key in seq.description for key in keys):
                    print(seq.format("fasta"))
                    to_write.append(seq)

        elif args.invert is True:
            if args.method == "exact":
                if seq.id not in keys:
                    print(seq.format("fasta"))
                    to_write.append(seq)
            elif args.method == "partial":
                if any(key not in seq.description for key in keys):
                    print(seq.format("fasta"))
                    to_write.append("fasta")

    if args.outfile:
        for rec in to_write:
            SeqIO.write(rec, args.outfile, "fasta")

if __name__ == "__main__":
    main()
