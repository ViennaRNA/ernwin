#!/usr/bin/python

from __future__ import absolute_import
from __future__ import print_function
import sys
from borgy.graph.bulge_graph import BulgeGraph
from borgy.utilities.vector import vec_distance
from borgy.graph.graph_pdb import get_bulge_centroid, get_mids

from numpy import dot

def get_c1_distance(chain, res1, res2):
    #print >>sys.stderr, "res1:", res1, "res2:", res2
    a1 = chain[res1]['P']
    a2 = chain[res2]['P']
    
    return a2 - a1
   
def print_bulge_distances(bg):
    for key in bg.edges:
        if key[0] != 's':
            edges = bg.edges[key]

            if len(bg.defines[key]) == 0:
                continue

            if len(bg.edges[key]) == 1:
                loop_length = abs(int(bg.defines[key][0]) - int(bg.defines[key][1]))
                distance = vec_distance(bg.coords[key][0], bg.coords[key][1])

                print("loop %d 0 %f" % ( loop_length, distance))

            elif len(bg.edges[key]) == 2:
                define = bg.defines[key]
                if len(define) == 4:
                    lengths = [abs(int(define[0]) - int(define[1])) + 1, abs(int(define[2]) - int(define[3])) + 1]
                elif len(define) == 2:
                    lengths = [abs(int(define[0]) - int(define[1])) + 1, 0]

                shortest_dist = vec_distance(bg.coords[key][0], bg.coords[key][1])
                lengths.sort()
                print("bulge_mid %d %d %f" % (lengths[0], lengths[1], shortest_dist))
                
def main():
    if len(sys.argv) < 2:
        print("Usage: ./graph_to_distances.py temp.comp")
        print()
        print("Traverse the bulge graph and print statistics on the distances between atoms")
        sys.exit(1)

    bg = BulgeGraph(sys.argv[1])
    print_bulge_distances(bg)

if __name__=="__main__":
    main()
