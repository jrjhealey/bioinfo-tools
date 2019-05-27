# A basic, and pretty fast (for python) way of generating all kmers of length x (only argument)
"""
Some speed benchmarks:
$ for i in {1..10..2} ; do
    echo "k =" $i
    time python mer-permutations.py $i > /dev/null
  done

k = 1
real    0m0.022s
user    0m0.016s
sys    0m0.004s

k = 3
real    0m0.022s
user    0m0.004s
sys    0m0.016s

k = 5
real    0m0.024s
user    0m0.016s
sys    0m0.004s

k = 7
real    0m0.054s
user    0m0.048s
sys    0m0.004s

k = 9
real    0m0.519s
user    0m0.504s
sys    0m0.016s
"""

import sys
import itertools

combinations = itertools.product(
    *itertools.repeat(["A", "T", "C", "G"], int(sys.argv[1]))
)
for i, k in enumerate(combinations):
    #    print('>Kmer_' + str(i) + '\n' + ''.join(k) )
    print("".join(k))
