# This script will calculate Shannon entropy and other
# information from a MSA.

# Dependencies:

# Biopython, Matplotlib, Math


#TODO
# 1. Add a brute force consensus sequence option (force choice between amino acids if equal likelihood (randomise?)
#    so that there are no ambiguous amino acids/bases
# 2. Do logo/conservation code (report back proportions of the MSA that correspond to a given base/amino).
# 3. Figure out the alphabet error for the information content function in Biopython.
# 4. Implement a mutation caller based on the following code:
	# import sys
	# from Bio import SeqIO
	#
	# f = sys.argv[1]
	#
	# seq_records = SeqIO.parse(f, 'fasta')
	# refseq_record = next(seq_records)
	# for seq_record in seq_records:
	#    for i in range(0, len(refseq_record)):
	#        nt1 = refseq_record[i]
	#        nt2 = seq_record[i]
	#        if nt1 != nt2:
	#            printseq_record.id, i+1, nt2, nt1)



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
__title__ = "MSAnalysis.py"
__license__ = "GPLv3"
__author_email__ = "J.R.J.Healey@warwick.ac.uk"


def parseArgs():
    """Parse command line arguments"""

    import argparse

    try:
        parser = argparse.ArgumentParser(
            description="Perform simple tasks with multiple sequence alignments \
			(e.g. Shannon Entropy, Consensus sequences, residue proportions). \
			Choose a task from [entropy|consensus|logo] and add any optional \
			parameters you want.")

	# Script task choices
	parser.add_argument('entropy',
			    help='Perform calculation of the Shannon Entropy (per column).')
	parser.add_argument('consensus',
			    nargs='+',
			    default=[0.51, ''],
			    help='Calculate a rough consensus sequence (from Bio.AlignIO). \'consensus\' should be followed by the threshold and ambiguous characters if needed (e.g. "consensus 0.5 N".')
	parser.add_argument('logo',
			    help='Calculate proportions of residues in each position of an MSA. Outputs a table to plot as a stacked bar/sequence logo.')
	# Required options
        parser.add_argument('-a',
                            '--alignment',
                            action='store',
                            required=True,
                            help='The multiple sequence alignment (MSA) in any of the formats supported by Biopython\'s AlignIO.')
        
        # Optional parameters
        parser.add_argument('-f',
                            '--alnformat',
                            action='store',
                            default='fasta',
                            help='Specify the format of the input MSA to be passed in to AlignIO.')
        parser.add_argument('-v',
                            '--verbose',
                            action='count',
                            default=1,
                            help='Verbose behaviour, printing parameters of the script and status messages.')
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



def parseMSA(msa, alnformat, consensus):
    """Parse in the MSA file using Biopython's AlignIO"""

    from Bio import AlignIO
    
    alignment = AlignIO.read(msa, alnformat)
    	
    # Do a little sanity checking:
    seq_lengths_list = []
    for record in alignment:
       seq_lengths_list.append(len(record))

    seq_lengths = set(seq_lengths_list)

    vprint("Alignment length is:" + str(list(seq_lengths)))

    if len(seq_lengths) != 1:
        sys.stderr.write("Your alignment lengths aren't equal. Check your alignment file.")
        sys.exit(1)

    index = range(1, list(seq_lengths)[0]+1)

    if consensus:
	from Bio.Align import AlignInfo
    	summary = AlignInfo.SummaryInfo(alignment)

    	vprint("Calculating a dumb consensus sequence with a threshold of: " + str(consensus[0] + " and taking " + str(consensus[1]) + " as the ambiguous character.")
    		
        consensus_seq = summary.dumb_consensus(consensus[0],consensus[1])



    return alignment, list(seq_lengths), index, consensus_seq

def getConsensus(alignment):
    """Uses summary functions from AlignIO to create a consensus sequence and perform operations on the MSA."""
    
    from Bio.Align import AlignInfo
    summary = AlignInfo.SummaryInfo(alignment)

    vprint("Calculating a dumb consensus sequence with a threshold of: "
            + str(consensus[0]
            + " and taking "
            + str(consensus[1])
            + " as the ambiguous character.")

    if consensus[1]:
	consensus_seq = summary.dumb_consensus(consensus[0],consensus[1])
    else:
	consensus_seq = summary.dumb_consensus(consensus[0])	

#info_content = summary.information_content()
    
    return consensus_seq
        
def shannon_entropy(list_input):
    """Calculate Shannon's Entropy per column of the alignment (H=-\sum_{i=1}^{M} P_i\,log_2\,P_i)"""

##################################################################
# Function to calcuate the Shannon's entropy per alignment column
# H=-\sum_{i=1}^{M} P_i\,log_2\,P_i (http://imed.med.ucm.es/Tools/svs_help.html)
# Gaps and N's are included in the calculation
##################################################################


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
    """Iterate across the whole MSA to calculate Shannon Entropy."""

    shannon_entropy_list = []
    for col_no in xrange(len(list(alignment[0]))):
        list_input = list(alignment[:, col_no])
        shannon_entropy_list.append(shannon_entropy(list_input))

    return shannon_entropy_list

def baseProportions(alignment):
    """Get the base proportions for each column of the alignment."""


def avgPairwise(alignment):
    """Calculate the average pairwise identity between sequences within an MSA. """    


def plot(index, sel):
    """"Create a quick plot via matplotlib to visualise the extended spectrum"""
    import matplotlib.pyplot as plt

    vprint("Plotting data...")

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


def vprint(verbose, string):
    if verbose > 0:
    	print(string)


#### MAIN CODE BEGINS ###
def main():
    """Compute Shannon Entropy from a provided MSA."""

    # Parse arguments
    args = parseArgs()

    # Convert object elements to standard variables for functions
    msa = args.alignment
    alnformat = args.alnformat
    global verbose
    verbose = args.verbose
    makeplot = args.makeplot
    runningmean = args.runningmean
    consensus = args.consensus
    logo = args.logo
    entropy = args.entropy

# Start calling functions to do the heavy lifting

    alignment, seq_lengths, index, consensus_seq = parseMSA(msa, alnformat, consensus, verbose)
    
    if entropy:
    	vprint("Performing entropy calculations on the provided MSA.")
    	sel = shannon_entropy_list_msa(alignment)

        if runningmean > 0:
            sel = running_mean(sel, runningmean)

	vprint("Index" + '\t' + "Entropy")
       	for c1, c2 in zip(index, sel):
       	    print(str(c1) + '\t' + str(c2))

	if makeplot is True:
	    plot(index, sel)

    elif consensus:
    	consensus_seq = getConsensus(alignment)
	from Bio import SeqIO
	SeqIO.write(consensus_seq, './consensus_sequence.fasta', fasta')    


    if runningmean > 0:
        sel = running_mean(sel, runningmean)

    if makeplot is True:
        plot(index, sel, verbose)


if __name__ == '__main__':
    main()

