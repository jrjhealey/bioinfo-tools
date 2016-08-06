cat effect.clust | grep Aligned | cut -d ' ' -f 6 | awk '{ sum +=  } END { print sum / NR }'
