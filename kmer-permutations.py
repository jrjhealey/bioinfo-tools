# A basic, and pretty fast (for python) way of generating all kmers of length x (only argument)
import sys
import itertools

combinations = itertools.product(*itertools.repeat(['A','T','C','G'], int(sys.argv[1])))
for i, k in enumerate(combinations):
    print('>Kmer_' + str(i) + '\n' + ''.join(k) )
