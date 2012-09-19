#!/usr/bin/python

import sys

from optparse import OptionParser
import pickle

from corgy.graph.bulge_graph import BulgeGraph
from corgy.builder.energy import LongRangeInteractionCount, CombinedEnergy
from corgy.builder.energy import DistanceIterator

from corgy.builder.models import SpatialModel

def main():
    parser = OptionParser()

    parser.add_option('-e', '--energy', dest='energy', default='', help="The energy function to use when evaluating structures")

    (options, args) = parser.parse_args()

    if len(args) < 1:
        print >>sys.stderr, "Usage: ./eval_energy.py temp.comp"
        print >>sys.stderr, "Evaluate the energy of a coarse-grained model"
        sys.exit(1)

    bg = BulgeGraph(args[0])
    bg.calc_bp_distances()
    sm = SpatialModel(bg)
    sm.get_sampled_bulges()

    default_energy = '/home/mescalin/pkerp/projects/ernwin/bobbins/energy/%s/1000/SkewNormalInteractionEnergy/LongRangeInteractionCount/JunctionClosureEnergy/CombinedEnergy.energy' % (bg.name)

    if options.energy == '':
        energy_function = pickle.load(open(default_energy, 'r'))
    else:
        energy_function = pickle.load(open(options.energy, 'r'))

    print energy_function.eval_energy(sm, True)


if __name__ == '__main__':
    main()

