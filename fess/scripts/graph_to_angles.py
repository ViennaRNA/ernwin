#!/usr/bin/python

import sys
from numpy import array, cross, dot
from math import pi
from Bio.PDB import calc_dihedral, Vector

from corgy.graph.bulge_graph import BulgeGraph
from corgy.utilities.vector import vec_distance, vec_angle, magnitude, vector_rejection
from corgy.utilities.vector import change_basis, get_standard_basis
from corgy.graph.graph_pdb import get_stem_orientation_parameters

def print_orientation_angles(bulge, connections):
    '''
    Print all of the angles that relate one stem to another.
    The angles needed are enumerated below:

    The vector actors are as follows:

    vec1[0] -> vec1[1] -> (bulge) -> vec2[0] -> vec2[1]
    twist1[0] -> twist1[1] -> (bulge) twist2[0] -> twist2[1]

    vec1_axis = vec1[1] - vec1[0]
    vec2_axis = vec2[1] - vec2[0]

    Angle 1: The angle to orient stem2 along the plane defined by the axis
        of stem1 and its near twist (twist[1])
        
        rotation around: mids1[1] - mids1[0]
        angle to be rotated: angle between twist1 and the rejection of vec2_axis
                             onto vec1_axis

    Angle 2: The angle to orient stem2 from position 1 to the same orientation
        as vector1

        rotation around: cross(vec1_axis, twist[1])
        angle: vec_angle(vec1_axis, transformed vec2_axis)

    Angle 3: The rotation needed to align twist1[1] with twist2[0]

        axis: vec1_axis
        angle: ?????

    @param bulge: The name of the bule region.
    @param connections: The names of the two stems that are connected by this bulge.
    '''

    s1 = connections[0]
    s2 = connections[1]

    s1d = bg.defines[s1]
    s2d = bg.defines[s2]
    bd = bg.defines[bulge]

    mids1 = bg.coords[s1]
    twists1 = bg.twists[s1]

    mids2 = bg.coords[s2]
    twists2 = bg.twists[s2]

    (s1b, s1e) = bg.get_sides(s1, bulge)
    (s2b, s2e) = bg.get_sides(s2, bulge)

    stem1_vec = mids1[s1b] - mids1[s1e]
    bulge_vec = mids2[s2b] - mids1[s1b]
    stem2_vec = mids2[s2e] - mids2[s2b]

    twists1_vec = [twists1[s1b], twists1[s1e]]
    twists2_vec = [twists2[s2e], twists2[s2b]]

    # Get the orientations for orienting these two stems
    (r, u, v) = get_stem_orientation_parameters(stem1_vec, twists_1_vec[0], stem2_vec)

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
