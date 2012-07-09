#!/usr/bin/python

from optparse import OptionParser
from bobbins_config import ConstructionConfig

from corgy.graph.bulge_graph import BulgeGraph
from corgy.utilities.data_structures import DefaultDict
from corgy.builder.models import SpatialModel
from corgy.builder.energy import LongRangeDistanceEnergy
from corgy.utilities.vector import get_vector_centroid, center_on_centroid
from corgy.builder.rmsd import rmsd, optimal_superposition

from corgy.builder.stats import AngleStatsDict, StemStatsDict

import copy
import sys
from sys import stderr

def main():
    if len(sys.argv) < 2:
        print >>sys.stderr, "Usage: ./construct.py temp.graph [constraint_file]"
        print
        print >>sys.stderr, "Construct an RNA molecule based on the graph"
        sys.exit(1)


    parser = OptionParser()

    parser.add_option('-a', '--angle_stats', dest='angle_stats_fn', default=ConstructionConfig.angle_stats_file, help='Location of the angle statistics file.') 
    parser.add_option('-i', '--iterations', dest='iterations', default=1, help='Number of structures to generate', type='int')
    parser.add_option('-t', '--output_type', dest='output_type', default='coord', help="Format of the output")
    parser.add_option('-b', '--best_filename', dest='best_filename', default='best.coord', help="The filename to dump the best (least rmsd structure) into", type='str')
    parser.add_option('-s', '--rot_segment', dest='rot_segment', default='s1', help="The name of segment around which to center all of the generated structures", type='str')
    parser.add_option('-x', '--text', dest='print_text', default=False, action='store_true', help='Print the names of the segments in the pymol output')

    (options, args) = parser.parse_args()

    angle_stats = AngleStatsDict(options.angle_stats_fn)
    stem_stats = StemStatsDict(options.angle_stats_fn)
    loop_stats = LoopStatsDict(options.angle_stats_fn)

    bg = BulgeGraph(args[0])
    orig_bg = copy.deepcopy(bg)

    sm = SpatialModel(bg, angle_stats, stem_stats, loop_stats)
    sm.rot_segment = options.rot_segment

    best_rmsd = 100000.0

    centers_orig = orig_bg.get_centers()
    centroid1 = get_vector_centroid(centers_orig)
    crds1 = center_on_centroid(centers_orig)
    orig_bg.translate_coords(-centroid1)
    orig_bg.output('orig.coords')

    energy_function = LongRangeDistanceEnergy()
    lowest_energy = 100000.0

    for i in range(options.iterations):
        sm.traverse_and_build()

        if options.output_type == 'pymol':
            sm.bg.output('this.coords')
            sm.pymol_printer.coordinates_to_pymol(sm.bg)
            #sm.pymol_printer.transform_segments(sm.first_translation, sm.first_rotation)
            sm.pymol_printer.transform_segments(sm.get_transform(options.rot_segment), sm.get_rotation(options.rot_segment))


        else:
            centers_new = sm.bg.get_centers()
            centroid2 = get_vector_centroid(centers_new)
            crds2 = center_on_centroid(centers_new)
            
            r = rmsd(crds1, crds2)
            
            if r < best_rmsd:
                sm.bg.translate_coords(-centroid2)
                rot = optimal_superposition(crds1, crds2)
                sm.bg.rotate_coords(rot)
                sm.bg.output(options.best_filename)
                best_rmsd = r

            print "rmsd:", r

    print >>stderr, "best_rmsd:", best_rmsd

    if options.output_type == 'pymol':
        sm.pymol_printer.print_text = options.print_text
        sm.pymol_printer.output_pymol_file()


if __name__ == "__main__":
    main()
