#!/usr/bin/python

import sys

from corgy.graph.bulge_graph import BulgeGraph
from corgy.utilities.vector import vec_distance

def output_long_range_distances(bg):
    for key1 in bg.longrange.keys():
        for key2 in bg.longrange[key1]:
            point1 = bg.get_point(key1)
            point2 = bg.get_point(key2)

            print "%s %s %f" % (key1, key2, vec_distance(point1, point2))

def main():
    if len(sys.argv) < 2:
        print >>sys.stderr, "Usage: ./long_range_distances.py temp.comp"
        sys.exit(1)

    bg = BulgeGraph(sys.argv[1])
    output_long_range_distances(bg)
    

if __name__ == '__main__':
    main()

