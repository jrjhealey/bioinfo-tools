cat effect.clust | grep Aligned | cut -d ' ' -f 6 | awk '{ sum += $1} END { print sum / NR }'
