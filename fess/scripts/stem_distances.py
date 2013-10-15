#!/usr/bin/python

import sys
from optparse import OptionParser

import itertools as it
import borgy.graph.bulge_graph as cgb
import borgy.utilities.vector as cuv
import borgy.graph.graph_pdb as cgg

def main():
    usage = './stem_distances.py temp.comp'
    usage += """
    List the distances between stems.
    """

    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    bg = cgb.BulgeGraph(args[0])

    for (s1, s2) in it.combinations(bg.stems(), 2):

        adjacent = 0
        if bg.are_any_adjacent_stems(s1, s2):
            adjacent = 1

        closest_points = cuv.line_segment_distance(bg.coords[s1][0],
                                                   bg.coords[s1][1],
                                                   bg.coords[s2][0],
                                                   bg.coords[s2][1])

        closest_distance = cuv.magnitude(closest_points[1] - closest_points[0])
        print "stem", closest_distance, adjacent, s1, s2

                        
if __name__ == '__main__':
    main()

