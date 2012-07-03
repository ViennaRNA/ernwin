#!/usr/bin/python

import sys

from corgy.graph.bulge_graph import BulgeGraph
from corgy.builder.energy import LongRangeInteractionCount
from corgy.builder.energy import DistanceIterator

def main():
    if len(sys.argv) < 2:
        print >>sys.stderr, "Usage: ./eval_energy.py temp.comp"
        print >>sys.stderr, "Evaluate the energy of a coarse-grained model"
        sys.exit(1)

    bg = BulgeGraph(sys.argv[1])
    lric_bounded = LongRangeInteractionCount()
    lric_naive = LongRangeInteractionCount(DistanceIterator())

    print lric_bounded.eval_energy(bg), lric_naive.eval_energy(bg)


if __name__ == '__main__':
    main()

