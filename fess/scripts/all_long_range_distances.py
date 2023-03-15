#!/usr/bin/python

from __future__ import absolute_import
from __future__ import print_function
import sys

from borgy.graph.bulge_graph import BulgeGraph
from borgy.utilities.vector import vec_distance

def output_all_distances(bg):
    for key1 in bg.defines.keys():
        for key2 in bg.defines.keys():
            if key1 != key2:
                if key2 not in bg.edges[key1]:
                    point1 = bg.get_point(key1)
                    point2 = bg.get_point(key2)

                    print("%s %s %f" % (key1, key2, vec_distance(point1, point2)))
def main():
    if len(sys.argv) < 2:
        print("Usage: ./long_range_distances.py temp.comp", file=sys.stderr)
        sys.exit(1)

    bg = BulgeGraph(sys.argv[1])
    output_all_distances(bg)
    

if __name__ == '__main__':
    main()

