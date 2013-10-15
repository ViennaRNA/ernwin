#!/usr/bin/python

import sys

from optparse import OptionParser
import pickle

import borgy.builder.energy as cbe
import borgy.builder.models as cbm
import borgy.graph.bulge_graph as cgb

def main():
    parser = OptionParser()

    parser.add_option('-e', '--energy', dest='energy', default='', help="The energy function to use when evaluating structures")

    (options, args) = parser.parse_args()

    if len(args) < 1:
        print >>sys.stderr, "Usage: ./eval_energy.py temp.comp"
        print >>sys.stderr, "Evaluate the energy of a coarse-grained model"
        sys.exit(1)

    energy = cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy()])
        
    bg = cgb.BulgeGraph(args[0])
    sm = cbm.SpatialModel(bg)
    sm.get_sampled_bulges()

    print energy.eval_energy(sm)


if __name__ == '__main__':
    main()

