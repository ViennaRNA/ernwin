#!/usr/bin/python

import sys

import pickle

from corgy.graph.bulge_graph import BulgeGraph
from corgy.builder.energy import LongRangeInteractionCount, CombinedEnergy
from corgy.builder.energy import DistanceIterator

from corgy.builder.models import SpatialModel

def main():
    if len(sys.argv) < 2:
        print >>sys.stderr, "Usage: ./eval_energy.py temp.comp"
        print >>sys.stderr, "Evaluate the energy of a coarse-grained model"
        sys.exit(1)

    bg = BulgeGraph(sys.argv[1])
    bg.calc_bp_distances()
    lric_bounded = LongRangeInteractionCount()
    lric_naive = LongRangeInteractionCount(DistanceIterator())
    jce = pickle.load(open('energies/jce.energy', 'r'))

    #print lric_bounded.eval_energy(bg), lric_naive.eval_energy(bg), jce.eval_energy(bg)


    ef1 = pickle.load(open('energies/%s/LongRangeInteractionCount' % (bg.name), 'r'))
    #ef2 = pickle.load(open('energies/%s/LongRangeDistanceEnergy' % (bg.name), 'r'))
    ef3 = pickle.load(open('energies/%s/JunctionClosureEnergy' % (bg.name), 'r'))
    ef4 = pickle.load(open('energies/%s/SkewNormalInteractionEnergy' % (bg.name), 'r'))

    #energy_function = CombinedEnergy([ef1, ef4, ef3])
    energy_function = CombinedEnergy([ef4])
    sm = SpatialModel(bg)
    sm.get_sampled_bulges()

    print energy_function.eval_energy(sm, True)


if __name__ == '__main__':
    main()

