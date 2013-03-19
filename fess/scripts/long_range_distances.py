#!/usr/bin/python

import sys
import itertools as it

from optparse import OptionParser

from corgy.graph.bulge_graph import BulgeGraph
from corgy.utilities.vector import vec_distance
import corgy.utilities.vector as cuv

from collections import defaultdict

def output_long_range_distances(bg):
    seen = defaultdict(set)

    for (key1, key2) in it.combinations(bg.longrange.keys()):
        (i1,i2) = cuv.line_segment_distance(bg.coords[key1][0], bg.coords[key1][1],
                                         bg.coords[key2][0], bg.coords[key2][1])

        #print "%s %s %f" % (key1, key2, vec_distance(point1, point2))
        seq1 = 'x'
        seq2 = 'x'

        if bg.get_type(key1) != 's' and bg.get_type(key1) != 'i':
            seq1 = bg.get_seq(key1)
        if bg.get_type(key2) != 's' and bg.get_type(key2) != 'i':
            seq2 = bg.get_seq(key2)

        print "%s %s %d %s %s %d %f %s %s %d" % (key1, 
                                     bg.get_type(key1), 
                                     bg.get_length(key1),
                                     key2, 
                                     bg.get_type(key2),
                                     bg.get_length(key2),
                                     cuv.magnitude(i2-i1),
                                     seq1, seq2, 1)

def output_all_distances(bg):
    for (key1, key2) in it.combinations(bg.defines.keys(), 2):
        longrange = 0

        if key2 in bg.longrange[key1]:
            longrange = 1

        #point1 = bg.get_point(key1)
        #point2 = bg.get_point(key2)

        (i1,i2) = cuv.line_segment_distance(bg.coords[key1][0], bg.coords[key1][1],
                                         bg.coords[key2][0], bg.coords[key2][1])

        #print "%s %s %f" % (key1, key2, vec_distance(point1, point2))
        seq1 = 'x'
        seq2 = 'x'

        if bg.get_type(key1) != 's' and bg.get_type(key1) != 'i':
            seq1 = bg.get_seq(key1)
        if bg.get_type(key2) != 's' and bg.get_type(key2) != 'i':
            seq2 = bg.get_seq(key2)

        print "%s %s %d %s %s %d %f %s %s %d" % (key1, 
                                     bg.get_type(key1), 
                                     bg.get_length(key1),
                                     key2, 
                                     bg.get_type(key2),
                                     bg.get_length(key2),
                                     cuv.magnitude(i2-i1),
                                     seq1, seq2, longrange)
def main():
    parser = OptionParser()

    parser.add_option('-a', '--all-pairs', dest='all_pairs', default=False, action='store_true', help='Print all interactions')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        print >>sys.stderr, "Usage: ./long_range_distances.py temp.comp"
        sys.exit(1)

    bg = BulgeGraph(args[0])

    if options.all_pairs:
        output_all_distances(bg)
    else:
        output_long_range_distances(bg)

    

if __name__ == '__main__':
    main()

