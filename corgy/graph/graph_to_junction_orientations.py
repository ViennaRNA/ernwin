#!/usr/bin/python

import sys
from bulge_to_fork_graph import parse_graph, BulgeGraph
from graph_to_simplified_pymol import get_bulge_centroid
from Bio.PDB import *
from helix_angle import get_mids, dot
from math import sqrt, pi
from shape_fit import vec_angle

def get_angle(chain, ds1, ds2):
    mids1 = get_mids(chain, ds1)
    mids2 = get_mids(chain, ds2)
    
    angle = vec_angle(mids1[1] - mids1[0], mids2[1] - mids2[0])
    return min(angle, pi - angle)

def print_helix_orientations(chain, bg):
    for d in bg.defines.keys():
        if d[0] == 's':
            continue

        if bg.weights[d] == 2:
            continue

        edges = bg.edges[d]
        for i in range(len(edges)):
            for j in range(i+1, len(edges)):
                if list(edges)[i][0] == 's' and list(edges)[j][0] == 's':
                    angle = get_angle(chain, bg.defines[list(edges)[i]], bg.defines[list(edges)[j]])
                    print "angle", bg.weights[d], angle
                
def main():
    if len(sys.argv) < 3:
        print "Usage: ./graph_to_distances.py temp.graph temp.pdb"
        print
        print "Traverse the bulge graph and print statistics on the distances between atoms"
        sys.exit(1)

    if sys.argv[1] == '-':
        f = sys.stdin
    else:
        f = open(sys.argv[1], 'r')

    chain = list(PDBParser().get_structure('temp', sys.argv[2]).get_chains())[0]

    bg = parse_graph(f)

    print_helix_orientations(chain, bg)

if __name__=="__main__":
    main()
