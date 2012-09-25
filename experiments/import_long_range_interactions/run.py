#!/usr/bin/python

import sys, copy
from optparse import OptionParser

import corgy.graph.bulge_graph as cgb
import corgy.builder.models as cbm
import corgy.builder.energy as cbe
import corgy.builder.sampling as cbs

def main():
    parser = OptionParser()

    parser.add_option('-i', '--iterations', dest='iterations', default=10, help='Number of structures to generate', type='int')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        print >>sys.stderr, "Usage: ./run.sh temp.comp"
        sys.exit(1)

    bg = cgb.BulgeGraph(args[0])
    bg.calc_bp_distances()
    sm = cbm.SpatialModel(bg)

    # Get the native long range interactions
    constraints = bg.get_long_range_constraints()

    for i in xrange(len(constraints)):
        # energy function consisting of just one native interaction
        energy_function = cbe.DistanceEnergy(constraints[i:i+1])
        plotter = None

        # keep track of statistics
        stats = cbs.SamplingStatistics(sm, plotter, 'b', silent=True)
        random_stats = cbs.SamplingStatistics(sm, plotter, 'r', silent=True)

        gs = cbs.GibbsBGSampler(copy.deepcopy(sm), energy_function, stats)
        gs_random = cbs.GibbsBGSampler(copy.deepcopy(sm), cbe.CombinedEnergy([cbe.RandomEnergy()]), random_stats)

        for j in range(options.iterations):
            gs.step()
            gs_random.step()

        sorted_energies = sorted(stats.energy_rmsd_structs, key=lambda x: x[0])

        node1 = constraints[i:i+1][0][0]
        node2 = constraints[i:i+1][0][1]

        print node1, node2, bg.bp_distances[node1][node2], sorted_energies[0][1]

if __name__ == '__main__':
    main()

