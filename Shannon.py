# This script will calculate Shannon entropy from a MSA.

# Dependencies:

# Biopython, Seaborn, Matplotlib, Math

"""
Shannon's entropy equation (latex format):
    H=-\sum_{i=1}^{M} P_i\,log_2\,P_i

    Entropy is a measure of the uncertainty of a probability distribution (p1, ..... , pM)
    https://stepic.org/lesson/Scoring-Motifs-157/step/7?course=Bioinformatics-Algorithms&unit=436

    Where, Pi is the fraction of nuleotide bases of nuleotide base type i,
    and M is the number of nuleotide base types (A, T, G or C)

    H ranges from 0 (only one base/residue in present at that position) to 4.322 (all 20 residues are equally
    represented in that position).

    Typically, positions with H >2.0 are considerered variable, whereas those with H < 2 are consider conserved.
    Highly conserved positions are those with H <1.0 (Litwin and Jores, 1992).
    A minimum number of sequences is however required (~100) for H to describe the diversity of a protein family.

"""
import os
import sys
import warnings
import traceback

__author__ = "Joe R. J. Healey"
__version__ = "1.0.0"
__title__ = "ShannonMSA"
__license__ = "GPLv3"
__author_email__ = "J.R.J.Healey@warwick.ac.uk"


def parseArgs():
    """Parse command line arguments"""

    import argparse

    try:
        parser = argparse.ArgumentParser(
            description='Compute per base/residue Shannon entropy of a Multiple Sequence Alignment.')

        parser.add_argument('-a',
                            '--alignment',
                            action='store',
                            required=True,
                            help='The multiple sequence alignment (MSA) in any of the formats supported by Biopython\'s AlignIO.')
        parser.add_argument('-f',
                            '--alnformat',
                            action='store',
                            default='fasta',
                            help='Specify the format of the input MSA to be passed in to AlignIO.')
        parser.add_argument('-v',
                            '--verbose',
                            action='count',
                            default=0,
                            help='Verbose behaviour, printing parameters of the script.')
        parser.add_argument('-m',
                            '--runningmean',
                            action='store',
                            type=int,
                            default=0,
                            help='Return the running mean (a.k.a moving average) of the MSAs Shannon Entropy. Makes for slightly smoother plots. Providing the number of points to average over switches this on.')
        parser.add_argument('--makeplot',
                            action='store_true',
                            help='Plot the results via Matplotlib.')
    except:
        print "An exception occurred with argument parsing. Check your provided options."
        traceback.print_exc()

    return parser.parse_args()



def parseMSA(msa, alnformat, verbose):
    """Parse in the MSA file using Biopython's AlignIO"""

    from Bio import AlignIO
    alignment = AlignIO.read(msa, alnformat)

    # Do a little sanity checking:
    seq_lengths_list = []
    for record in alignment:
       seq_lengths_list.append(len(record))

    seq_lengths = set(seq_lengths_list)

    if verbose > 0: print("Alignment length is:" + str(list(seq_lengths)))

    if len(seq_lengths) != 1:
        sys.stderr.write("Your alignment lengths aren't equal. Check your alignment file.")
        sys.exit(1)

    index = range(1, list(seq_lengths)[0]+1)

    return alignment, list(seq_lengths), index

##################################################################
# Function to calcuate the Shannon's entropy per alignment column
# H=-\sum_{i=1}^{M} P_i\,log_2\,P_i (http://imed.med.ucm.es/Tools/svs_help.html)
# Gaps and N's are included in the calculation
##################################################################

def shannon_entropy(list_input):
    """Calculate Shannon's Entropy per column of the alignment (H=-\sum_{i=1}^{M} P_i\,log_2\,P_i)"""

    import math
    unique_base = set(list_input)
    M   =  len(list_input)
    entropy_list = []
    # Number of residues in column
    for base in unique_base:
        n_i = list_input.count(base) # Number of residues of type i
        P_i = n_i/float(M) # n_i(Number of residues of type i) / M(Number of residues in column)
        entropy_i = P_i*(math.log(P_i,2))
        entropy_list.append(entropy_i)

    sh_entropy = -(sum(entropy_list))

    return sh_entropy


def shannon_entropy_list_msa(alignment):
    """Calculate Shannon Entropy across the whole MSA"""

    shannon_entropy_list = []
    for col_no in xrange(len(list(alignment[0]))):
        list_input = list(alignment[:, col_no])
        shannon_entropy_list.append(shannon_entropy(list_input))

    return shannon_entropy_list


def plot(index, sel, verbose):
    """"Create a quick plot via matplotlib to visualise the extended spectrum"""
    import matplotlib.pyplot as plt

    if verbose > 0: print("Plotting data...")

    plt.plot(index, sel)
    plt.xlabel('MSA Position Index', fontsize=16)
    plt.ylabel('Shannon Entropy', fontsize=16)

    plt.show()


def running_mean(l, N):
    sum = 0
    result = list(0 for x in l)

    for i in range( 0, N ):
        sum = sum + l[i]
        result[i] = sum / (i+1)

    for i in range( N, len(l) ):
        sum = sum - l[i-N] + l[i]
        result[i] = sum / N

    return result

def main():
    """Compute Shannon Entropy from a provided MSA."""

    # Parse arguments
    args = parseArgs()

    # Convert object elements to standard variables for functions
    msa = args.alignment
    alnformat = args.alnformat
    verbose = args.verbose
    makeplot = args.makeplot
    runningmean = args.runningmean

# Start calling functions to do the heavy lifting

    alignment, seq_lengths, index = parseMSA(msa, alnformat, verbose)
    sel = shannon_entropy_list_msa(alignment)

    if runningmean > 0:
        sel = running_mean(sel, runningmean)

    if makeplot is True:
        plot(index, sel, verbose)

    if verbose > 0: print("Index" + '\t' + "Entropy")
    for c1, c2 in zip(index, sel):

        print(str(c1) + '\t' + str(c2))



if __name__ == '__main__':
    main()

