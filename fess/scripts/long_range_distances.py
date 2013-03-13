#!/usr/bin/python

import sys

from corgy.graph.bulge_graph import BulgeGraph
from corgy.utilities.vector import vec_distance
import corgy.utilities.vector as cuv

from collections import defaultdict

def output_long_range_distances(bg):
    seen = defaultdict(set)

    for key1 in bg.longrange.keys():
        for key2 in bg.longrange[key1]:
            if key2 in seen[key1]:
                continue

            seen[key1].add(key2)
            seen[key2].add(key1)

            #point1 = bg.get_point(key1)
            #point2 = bg.get_point(key2)

            (i1,i2) = cuv.line_segment_distance(bg.coords[key1][0], bg.coords[key1][1],
                                             bg.coords[key2][0], bg.coords[key2][1])

            #print "%s %s %f" % (key1, key2, vec_distance(point1, point2))
            print "%s %s %f" % (key1, key2, cuv.magnitude(i2-i1))

def main():
    if len(sys.argv) < 2:
        print >>sys.stderr, "Usage: ./long_range_distances.py temp.comp"
        sys.exit(1)

    bg = BulgeGraph(sys.argv[1])
    output_long_range_distances(bg)
    

if __name__ == '__main__':
    main()

