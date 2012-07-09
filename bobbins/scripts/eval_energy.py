#!/usr/bin/python

import sys

import pickle

from corgy.graph.bulge_graph import BulgeGraph
from corgy.builder.energy import LongRangeInteractionCount, CombinedEnergy
from corgy.builder.energy import DistanceIterator

def main():
    if len(sys.argv) < 2:
        print >>sys.stderr, "Usage: ./eval_energy.py temp.comp"
        print >>sys.stderr, "Evaluate the energy of a coarse-grained model"
        sys.exit(1)

    bg = BulgeGraph(sys.argv[1])
    lric_bounded = LongRangeInteractionCount()
    lric_naive = LongRangeInteractionCount(DistanceIterator())

    '''
    ef1 = pickle.load(open('energies/lric.energy', 'r'))
    ef2 = pickle.load(open('energies/lrde.energy', 'r'))

    energy_function = CombinedEnergy([ef1, ef2])
    '''

    print lric_bounded.eval_energy(bg), lric_naive.eval_energy(bg)
    #print energy_function.eval_energy(bg)


if __name__ == '__main__':
    main()

