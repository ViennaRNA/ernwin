#!/usr/bin/python

import sys
from optparse import OptionParser

import corgy.graph.bulge_graph as cgb
import corgy.graph.graph_pdb as cgg

import itertools as it

def stem_stem_orientations(bg):
    for (s1, s2) in it.combinations(bg.stems(), r=2):
        if not bg.are_adjacent_stems(s1, s2):
            print " ".join(map(str, cgg.stem_stem_orientation(bg, s1, s2)))

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

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    bg = cgb.BulgeGraph(args[0])
    stem_stem_orientations(bg)

if __name__ == '__main__':
    main()

