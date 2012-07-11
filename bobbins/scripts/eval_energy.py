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
    jce = pickle.load(open('energies/jce.energy', 'r'))

    #print lric_bounded.eval_energy(bg), lric_naive.eval_energy(bg), jce.eval_energy(bg)

    ef1 = pickle.load(open('energies/lric.energy', 'r'))
    ef2 = pickle.load(open('energies/lrde.energy', 'r'))
    ef3 = pickle.load(open('energies/jce.energy', 'r'))
    ef4 = pickle.load(open('energies/SkewNormalInteractionEnergy.energy', 'r'))

    energy_function = CombinedEnergy([ef1, ef4, ef3])


    print energy_function.eval_energy(bg, True)


if __name__ == '__main__':
    main()

