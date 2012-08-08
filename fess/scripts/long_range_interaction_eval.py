#!/usr/bin/python

import sys


from corgy.graph.bulge_graph import BulgeGraph
from corgy.builder.energy import LongRangeInteractionCount, CombinedEnergy
from corgy.builder.energy import DistanceIterator

def main():
    if len(sys.argv) < 2:
        print >>sys.stderr, "Usage: ./long_range_interaction_count.py temp.comp"
        sys.exit(1)

    bg = BulgeGraph(sys.argv[1])
    lric_bounded = LongRangeInteractionCount()
    lric_naive = LongRangeInteractionCount(DistanceIterator())

    print lric_bounded.count_interactions(bg), lric_naive.count_interactions(bg)

if __name__ == '__main__':
    main()

