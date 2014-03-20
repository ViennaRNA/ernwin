#!/usr/bin/python

import sys
from optparse import OptionParser

import forgi.threedee.model.coarse_grain as ttmc
import borgy.graph.graph_pdb as cgg

import itertools as it

def stem_stem_orientations(bg):
    for (s1, s2) in it.permutations(bg.stem_iterator(), r=2):
        if not bg.are_adjacent_stems(s1, s2):
            print " ".join(map(str, cgg.stem_stem_orientation(bg, s1, s2)))

def loop_loop_orientations(bg):
    for (l1, l2) in it.permutations(bg.loops(), r=2):
        if l1 != l2:
            print " ".join(map(str, cgg.stem_stem_orientation(bg, l1, l2)))

def main():
    usage = './stem_stem_orientations.py temp.comp'
    usage += '''
    Output the orientations of the stems as a set of'
    three coordinates.
    1. The distance between the closest points of the
    two stems.
    2. The angle between the two stems in the plane defined
    by the axis of stem1 and the vector between the closest
    points between the two stem.
    3. The angle of the second stem out the plane.
    '''

    parser = OptionParser()

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')
    parser.add_option('-l', '--loops', dest='loops', default=False, action='store_true', help="Compute the statistics for the loop regions rather than the stems.")

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    for arg in args:
        bg = ttmc.CoarseGrainRNA(arg)
        if options.loops:
            loop_loop_orientations(bg)
        else:
            stem_stem_orientations(bg)

if __name__ == '__main__':
    main()

