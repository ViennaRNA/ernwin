#!/usr/bin/python

import sys
from numpy import array, cross, dot
from math import pi
from Bio.PDB import calc_dihedral, Vector

from corgy.graph.bulge_graph import BulgeGraph
from corgy.utilities.vector import vec_distance, vec_angle, magnitude, vector_rejection
from corgy.utilities.vector import change_basis, get_standard_basis
from corgy.graph.graph_pdb import get_stem_orientation_parameters
from corgy.graph.graph_pdb import get_stem_separation_parameters
from corgy.graph.graph_pdb import print_orientation_angles

from corgy.graph.graph_pdb import get_stem_twist_and_bulge_vecs
from corgy.graph.graph_pdb import get_twist_angle

def print_new_bulge_angles(bg):
    '''
    This should be refactored to use the corgy.builder.stats.AngleStat class.
    '''
    for define in bg.defines.keys():
        if define[0] != 's' and len(bg.edges[define]) == 2:
            (as1, as2) = bg.get_bulge_angle_stats(define)

            print "angle", as1.pdb_name, as1.dim1, as1.dim2, as1.u, as1.v, as1.t, as1.r1, as1.u1, as1.v1
            print "angle", as2.pdb_name, as2.dim1, as2.dim2, as2.u, as2.v, as2.t, as2.r1, as2.u1, as2.v1
        else:
            continue

def print_stem_stats(bg):
    for d in bg.defines.keys():
        if d[0] == 's':
            base_pair_length = abs(bg.defines[d][0] - bg.defines[d][1])
            phys_length = magnitude(bg.coords[d][1] - bg.coords[d][0])
            twist_angle = get_twist_angle(bg.coords[d], bg.twists[d])

            print "stem", bg.name, base_pair_length, phys_length, twist_angle

def print_loop_stats(bg):
    for d in bg.defines.keys():
        if d[0] != 's':
            if bg.weights[d] == 1 and len(bg.edges[d]) == 1:
                # unpaired region
                base_pair_length = abs(bg.defines[d][0] - bg.defines[d][1])
                phys_length = magnitude(bg.coords[d][1] - bg.coords[d][0])
                
                print "loop", bg.name, base_pair_length, phys_length

def main():
    if len(sys.argv) < 2:
        print >>sys.stderr, "Usage: ./graph_to_angles.py temp.comp"
        print
        print >>sys.stderr, "Traverse a structure and output the stems that are connected by a bulge"
        run_tests()
        sys.exit(1)

    if sys.argv[1] == '-':
        f = sys.stdin
    else:
        f = open(sys.argv[1], 'r')

    bg = BulgeGraph(sys.argv[1])
    print_new_bulge_angles(bg)
    print_stem_stats(bg)
    print_loop_stats(bg)

if __name__=="__main__":
    main()
