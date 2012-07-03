#!/usr/bin/python

from optparse import OptionParser

from bobbins_config import ConstructionConfig

from corgy.graph.bulge_graph import BulgeGraph
from corgy.builder.stats import AngleStatsDict, StemStatsDict
from corgy.builder.energy import LongRangeDistanceEnergy
from corgy.builder.models import SpatialModel
from corgy.builder.rmsd import centered_rmsd

from corgy.utilities.vector import get_vector_centroid, center_on_centroid

import sys
from sys import stderr

def main():
    parser = OptionParser()

    parser.add_option('-e', '--energy', dest='energy', default='lrde', help="The energy function to use when evaluating structures")
    parser.add_option('-i', '--iterations', dest='iterations', default=10, help='Number of structures to generate', type='int')
    parser.add_option('-b', '--best_filename', dest='best_filename', default='best.coord', help="The filename to dump the best (least rmsd structure) into", type='str')
    parser.add_option('-a', '--angle_stats', dest='angle_stats_fn', default=ConstructionConfig.angle_stats_file, help='Location of the angle statistics file.') 

    (options, args) = parser.parse_args()

    if len(args) < 1:
        print "Usage: ./sample.py temp.comp"
        sys.exit(1)

    angle_stats = AngleStatsDict(options.angle_stats_fn)
    stem_stats = StemStatsDict(options.angle_stats_fn)

    bg = BulgeGraph(args[0])
    sm = SpatialModel(bg, angle_stats, stem_stats)

    '''
    Get the vectors of the elements of the original
    structure.
    '''
    centers_orig = bg.get_centers()
    energy_function = LongRangeDistanceEnergy()
    lowest_energy = 100000.0

    energy_function.calibrate(bg)

    # print the native energy
    energy = energy_function.eval_energy(bg)
    print >>stderr, "original native_energy:", energy

    for i in range(options.iterations):
        sm.traverse_and_build()
        energy = energy_function.eval_energy(sm.bg)

        centers_new = sm.bg.get_centers()
        r = centered_rmsd(centers_orig, centers_new)

        print "native_energy:", energy, r

        if energy < lowest_energy:
            lowest_energy = energy
            lowest_rmsd = r
            sm.bg.output(options.best_filename)

    print >>stderr, "lowest energy:", lowest_energy, lowest_rmsd

if __name__ == '__main__':

    main()

