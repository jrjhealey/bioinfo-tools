#! /usr/bin/python

"""
MergeAlign algorithm for create consensus multiple sequence alignments

For a lsit of options, run:
>> python mergeAlign.py -h
"""

import os
import sys
import getopt

class Node:
    """ A position in alignment space with edges to the previous nodes in all alignments """
    
    def __init__(self):
        self.previous_nodes = {}
        self.path_score = 0
        self.path_length = 0
        self.path_average = 0
        self.best_previous_node = None

def combineMultipleAlignments(alignment_names):
    """ Combine all alignments in a folder into a single alignment.
        Takes a list filepaths to FASTA alignments.
        Returns a consensus alignment as a dictionary and a list of column scores. """
    
    # Read alignments
    alignments = [parseFASTA(alignment) for alignment in alignment_names]
    
    # Open the first file and remove gaps to get original sequences in a dictionary
    original_sequences = dict([(name, seq.replace('-', '')) for name, seq in alignments[0].iteritems()])
    sequence_names = original_sequences.keys()
    
    # Convert dictionaries of sequences to lists of indices
    indices = [[convertSequenceToIndices(alignment[seq]) for seq in sequence_names] for alignment in alignments]
    
    # Convert alignments into list of coordinates
    coordinates = [zip(*alignment) for alignment in indices]

    # Move through coordinates, generating paths of nodes through the alignments and finding best alignment
    nodes = createNodes(coordinates)
    final_coordinates, scores = scoreNodes(nodes, len(coordinates))
    final_alignment = convertListOfCoordinatesToSequences(final_coordinates, original_sequences)
    
    return final_alignment, scores

def parseFASTA(filename):
    """ Read file in FASTA format.
        Return a dictionary with key=sequence name, value=sequence. """

    try:
        f = open(filename, 'r')
    except IOError:
        print "Unable to open file %s" % filename
        sys.exit(2)

    sequences = {}
    seq_id = None

    for line in f.readlines():
        if line.startswith('>'):
            seq_id = line.rstrip()[1:]
            sequences[seq_id] = ''
        elif seq_id != None:
            sequences[seq_id] += line.rstrip()

    if len(sequences) == 0:
        print "%s contains no sequences" % filename
        sys.exit(2)

    return sequences

def convertSequenceToIndices(sequence):
    """ Convert a sequence to a list of amino acid indices.
        Gap are converted to the index of the preceeding amino acid. """
    
    indices = []
    i = 0
    
    for aa in sequence:
        if aa != '-':
            i += 1
        indices.append(i)
        
    return indices

def convertIndicesToSequence(indices, original_sequence):
    """ Convert list of indices back to original sequence plus gaps """
    
    sequence = ''

    previous_i = 0
    for i in indices:
        if i != previous_i:
            sequence += original_sequence[i-1]
        else:
            sequence += '-'
        previous_i = i
        
    return sequence

def convertListOfCoordinatesToSequences(final_coordinates, original_sequences):
    """ Convert a list of coordinates in alignment space into an alignment using a passed dictionary of sequences. """
    
    seq_names = original_sequences.keys()
    alignment = {}
    
    for i, sequence in enumerate(zip(*final_coordinates)):
        seq = seq_names[i]
        alignment[seq] = convertIndicesToSequence(sequence, original_sequences[seq])
    
    return alignment

def createNodes(alignments):
    """ Traverse each alignment creating nodes at each point found in alignment space. 
        For each node record which nodes preceded it and how often. """
    
    dimensions = len(alignments[0][0])
    first_node = tuple([0]*dimensions)
    nodes = {first_node: Node()}
    
    for alignment in alignments:
        previous_node = first_node
        
        for point in alignment:
            if not point in nodes:
                nodes[(point)] = Node()
            node_dict = nodes[(point)].previous_nodes
            # Add previous node to current nodes count of previous nodes
            node_dict[previous_node] = node_dict.get(previous_node, 0) + 1
            previous_node = point

    return nodes

def scoreNodes(nodes, num_paths=100):
    """ Travserse nodes, finding the best route to each one with a dynamic programming approach. 
        Return the best path to the final node. """
    
    dimensions = len(nodes.keys()[0])
    first_node = tuple([0]*dimensions)
    
    for coord, current_node in sorted(nodes.iteritems())[1:]:
        previous_nodes = current_node.previous_nodes.items()
        best_node = max(previous_nodes, key=lambda n: n[1] + nodes[n[0]].path_average)
        
        # Update this node with the best path so far
        current_node.path_length = nodes[best_node[0]].path_length + 1
        current_node.path_score = best_node[1] + nodes[best_node[0]].path_score
        current_node.path_average = current_node.path_score* 1.0 / current_node.path_length
        current_node.best_previous_node = best_node[0]
    
    #print "final score: %.2f; length: %d" % (current_node.path_average/num_paths, current_node.path_length)
    
    # Start at the final node and work backwards, moving from each node to its best previous node
    path = []
    scores = []
    while coord != first_node:
        path.append(coord)
        scores.append(nodes[coord].previous_nodes[nodes[coord].best_previous_node]*1.0/num_paths)
        coord = nodes[coord].best_previous_node

    path.reverse()
    scores.reverse()
    
    return path, scores

def outputAlignmentAsFASTA(filename, final_alignment, threshold=None, scores=None):
    """ Output alignment in form of a dictionary as a FASTA file. """
    
    output = file(filename, 'w')
    
    for name, sequence in final_alignment.items():
        output.write(">%s\n" % name)
        if threshold:
            sequence = ''.join([aa for (aa, score) in zip (sequence, scores) if score>threshold])
        output.write("%s\n" % sequence)
        
def outputAlignmentAsCLUSTAL(filename, final_alignment, scores):
    """ Output dictionary of alignments as a FASTA file. """
    
    output = file(filename, 'w')
    output.write("CLUSTAL format alignment by MultiMatrixMafft\n\n")
    
    trunc_names = [len(key)>19 and "%s..." % key[:16] or key for key in final_alignment.keys()]
    
    start = 0
    while start < len(scores):
        stop = start+60
        for i, sequence in enumerate(final_alignment.values()):
            output.write(">%s\t" % trunc_names[i])
            output.write("%s\n" % sequence[start:stop])
            
        output.write(" " * 24)        
        for score in scores[start:stop]:
            if score == 1:
                output.write("*")
            elif score > 0.9:
                output.write(":")
            elif score > 0.75:
                output.write(".")
            else:
                output.write(" ")
                
        output.write('\n\n')
        
        start = stop

def outputScore(filename, scores):
    """ Output score as a list of numbers. """
    
    output = file(filename, 'w')
    for score in scores:
        output.write('%.3f\n' % score)

def outputAlignment(filename, final_alignment, scores):    
    output = file(filename, 'w')
    
    for score in scores:
        if score == 1:
            output.write('*')
        else:
            output.write('%d' % (10*score))
    output.write('\n')
    
    for name, sequence in final_alignment.items():
        output.write(">%s\n" % name)
        output.write("%s\n" % sequence)

def outputAlignmentVertical(filename, final_alignment, original_sequences):
    output = file(filename, 'w')
    
    for i, sequence in enumerate(zip(*[residue for (residue, score) in final_alignment])):
        print convertIndicesToSequence(sequence, original_sequences[i])

def usage():
    print """
To combine multiple alignments:
    % python multimatrixmafft.py -a path/to/folder
    
  -a --alignments
    ARG: path to a folder
    Path to a folder containing all the FASTA alignments you want to combine.
    All the alignments must be of the same set of sequences.
    
  OPTIONAL 
  -f --fasta
    ARG: filename
    Final alignment is saved as a FASTA file with this name.
    
  -s --score
    ARG: filename
    Final scores are saved as a list in a text file with this name.
    
  -t --score
    ARG: number >0 and <= 1
    Only columns with a score > threshold are outputted.
    Only works in conjection with -f.

Additional arguments:

  -h, --help
    ARG: none
    Get command line arguments.
"""

if __name__ == '__main__':
    # Options
    matrices = None
    sequence_file = None
    alignment_folder = None
    fasta_output = None
    clustal_output = None
    score_output = None
    threshold = None
    
#   Get commandline arguments
    try:                                
        opts, args = getopt.getopt(sys.argv[1:],
                                   "ha:f:s:t:",
                                   ["help", "alignments=", "fasta=", "score=", "threshold="])
    except getopt.GetoptError:
        print "Error: command line argument not recognised"
        usage()
        sys.exit(2)
        
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        
        elif opt in ("-a", "--alignments"):
            alignment_folder = arg
        
        elif opt in ("-f", "--fasta"):
            fasta_output = arg
        
        elif opt in ("-s", "--score"):
            score_output = arg
            
        elif opt in ("-t", "--threshold"):   
            try:
                threshold = float(arg)
            except ValueError:
                print 'Error: threshold must be a number'
                usage()
                sys.exit(2)
            if not 0 < threshold <= 1:
                print 'Error: threshold must be >0 and <= 1'
                usage()
                sys.exit(2)
    
    # Ensure there is at least one argument
    if not alignment_folder:
        print "Error: alignment folder not defined"
        usage()
        sys.exit(2)

    # Get list of dictionaries of aligned sequences; indices correspond to matrices
    alignment_names = [os.path.join(alignment_folder, alignment) for alignment in os.listdir(alignment_folder)]
    # Combine alignments
    final_alignment, scores = combineMultipleAlignments(alignment_names)
    
    # Output alignment
    if fasta_output:
        outputAlignmentAsFASTA(fasta_output, final_alignment, threshold, scores)
    if score_output:
        outputScore(score_output, scores)