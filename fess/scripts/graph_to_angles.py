#!/usr/bin/python

import sys
from numpy import array, cross, dot
from math import pi
from Bio.PDB import calc_dihedral, Vector

from corgy.graph.bulge_graph import BulgeGraph
from corgy.utilities.vector import vec_distance, vec_angle, magnitude, vector_rejection
from corgy.utilities.vector import change_basis, get_standard_basis
from corgy.graph.graph_pdb import get_stem_orientation_parameters
from corgy.graph.graph_pdb import print_orientation_angles

def print_new_bulge_angles(bg):
    for define in bg.defines.keys():
        if define[0] != 's' and len(bg.edges[define]) == 2:
            connections = list(bg.edges[define])

            print_orientation_angles(define, connections)
        else:
            continue

def print_bulge_angles(bg):
    '''
    Print statistics of the angles.
    '''

    bulges = bg.get_bulge_vertices()

    for define in bg.defines.keys():
        # Only look at statistics for bulges which connect two stems

        if define[0] != 's' and len(bg.edges[define]) == 2:
            connections = list(bg.edges[define])
        else:
            continue

        bulge = define

        s1 = connections[0]
        s2 = connections[1]

        s1d = bg.defines[s1]
        s2d = bg.defines[s2]
        bd = bg.defines[bulge]

        # the starts and ends of the helices
        mids1 = bg.coords[s1]
        mids2 = bg.coords[s2]

        #get the side of each stem that is closest to the bulge
        (s1b, s1e) = bg.get_sides(s1, bulge)
        (s2b, s2e) = bg.get_sides(s2, bulge)
           
        out_str = "angle "

        if len(bd) == 2:
            # Bulges that are single stranded
            dims = (0, abs(bd[1] - bd[0]))
        else:
            # Double stranded bulges
            dims = (abs(bd[1] - bd[0]), abs(bd[2] - bd[3]))

        # s1vec -> bv -> s2vec
        stem1_vec = mids1[s1e] - mids1[s1b]
        bulge_vec = mids1[s1b] - mids2[s2b]
        stem2_vec = mids2[s2b] - mids2[s2e]

        # aliases for functions which require an array instead of a vector
        stem1_array = stem1_vec
        bulge_vec_array = bulge_vec

        s1b_cross = cross(bulge_vec_array, stem1_array)

        angle1 = vec_angle(stem1_vec, bulge_vec)
        angle2 = vec_angle(stem2_vec, bulge_vec)

        stem2_along_stem1 = dot(stem1_vec, stem2_vec) / magnitude(stem1_vec)
        stem2_along_s1b_cross = dot(s1b_cross, stem2_vec) / magnitude(s1b_cross)

        '''
        if stem2_along_s1b_cross < 0.:
            angle3 = vec_angle(array([1., 0., 0.]), array([stem2_along_stem1, stem2_along_s1b_cross, 0.]))
        else:
            angle3 = -vec_angle(array([1., 0., 0.]), array([stem2_along_stem1, stem2_along_s1b_cross, 0.]))
        '''
        angle3 = vec_angle(array([1., 0., 0.]), array([stem2_along_stem1, stem2_along_s1b_cross, 0.]))

        angle4 = vec_angle(stem1_vec, stem2_vec)

        angle5 = calc_dihedral(Vector(mids1[s1e]), Vector(mids1[s1b]), Vector(mids2[s2b]), Vector(mids2[s2e]))
        
        #angle3 = vec_angle(stem2_vec, s1b_cross)

        #print >>sys.stderr, "stem1_ncl:", stem1_ncl
        #print >>sys.stderr, "bulge_vec:", bulge_vec

        # m1(s1e) - m1(s1b) - m2(s2b) - m2(s2e)
        # s1vec -> bv -> s2vec

        out_str5 = "angles %s %d %d %f %f %f %f %f" % (bg.name, min(dims), max(dims), angle1, angle2, angle3, angle4, angle5)

        print out_str5

    #print bg.get_named_define

def run_tests():
    '''
    Run some tests for the angle caculation.
    '''

    vec1 = array([0., 5., 0.])
    vec2 = array([5., 5., 0.])

    twist1 = array([0., 0., 1.])

    print get_rotation_matrix_and_angles(array([0., 5., 0.]), array([5., 5., 0.]), array([0., 0., 1.]))
    print get_rotation_matrix_and_angles(array([0., 5., 0.]), array([5., 5., 0.]), array([0., 0., -1.]))

    print
    print

    print get_rotation_matrix_and_angles(array([0., 5., 0.]), array([5., 5., 0.]), array([0., -5., 1.]))
    print get_rotation_matrix_and_angles(array([0., 5., 0.]), array([5., 5., 0.]), array([0., -5., -1.]))

    print
    print

    print get_rotation_matrix_and_angles(array([2., 1., 1.]), array([0., 5., 0.]), array([0., 0., 1.]))
    print get_rotation_matrix_and_angles(array([2., 1., 1.]), array([0., 5., 0.]), array([0., 0., -1.]))

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

if __name__=="__main__":
    main()
