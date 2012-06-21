#!/usr/bin/python

import unittest
from sys import stderr

from random import random, uniform

from numpy import array, matrix, allclose, cross
from corgy.utilities.vector import change_basis, get_standard_basis, normalize
from corgy.utilities.vector import spherical_cartesian_to_polar, spherical_polar_to_cartesian
from corgy.utilities.vector import rotation_matrix
from corgy.utilities.vector import get_non_colinear_unit_vector
from math import pi, sqrt, asin

from Bio.PDB import PDBParser
from numpy import dot
from corgy.graph.graph_pdb import get_mids, get_twists
from corgy.graph.graph_pdb import get_stem_orientation_parameters
from corgy.graph.bulge_graph import BulgeGraph


class TestGraphToAngles(unittest.TestCase):
    '''
    Tests for the gathering of angles statistics.
    '''

    def orientation_test(self, stem1, twist1, stem2):
        print "\n\n"
        stem1 = normalize(array(stem1))
        stem2 = normalize(array(stem2))
        twist1 = normalize(array(twist1))

        (r, u, v) = get_stem_orientation_parameters(stem1, twist1, stem2)

        print "stem1:", stem1
        print "stem2:", stem2
        print "twist1:", stem2
        print "r:", r, "u:", u, "v:", v

        # rotate around the z-axis to align with the x-axis
        rot_mat1 = rotation_matrix(array([0., 0., 1.]), -v)

        # rotate around the y-axis to get back to [0., 0. 0.]
        rot_mat2 = rotation_matrix(array([0., 1., 0.]), u - pi/2) 

        stem2_new = dot(rot_mat1, stem2)
        print "stem2_new:", stem2_new
        stem2_new = dot(rot_mat2, stem2_new)
        print "stem2_new:", stem2_new

        self.assertTrue(allclose(stem2_new, stem1))
    


    def test_get_orientation_parameters(self):
        #stem1 = array([1., 0., 0.])
        #stem2 = array([1., 0., 0.])
        #twist1 = array([0., 1., 0.])

        #self.orientation_test(stem1, stem2, twist1)
        self.orientation_test([1., 0., 0.], [0., -1., 0.], [1, 0, 0])
        self.orientation_test([1., 0., 0.], [0., 1., 0.], [1, 0, 0])
        self.orientation_test([1., 1., 0.], [1., -1., 0.], [1, 0, 0])
        self.orientation_test([2., 1., 0.], [1., -1., 0.], [1, 0, 0])

        return

        for i in range(1):
                stem1 = [uniform(-1., 1.), uniform(-1., 1.), uniform(-1., 1.)]
                stem2 = [uniform(-1., 1.), uniform(-1., 1.), uniform(-1., 1.)]
                twist1 = cross(stem1, get_non_colinear_unit_vector(stem1))

                self.orientation_test(stem1, twist1, stem2)


class TestGraphPDBFunctions(unittest.TestCase):
    '''
    Tests for the functions in corgy.graph.graph_pdb.
    '''

    def test_mids_and_twists(self):
        '''
        The twist vectors should always be perpendicular to the stem vectors.
        '''

        s = PDBParser().get_structure('test', 'test/structures/1gid.pdb')
        chain = list(s.get_chains())[0]

        bg = BulgeGraph('test/graph/1gid.comp')

        for d in bg.defines.keys():
            if d[0] == 's':
                mids = get_mids(chain, bg.defines[d])
                twists = get_twists(chain, bg.defines[d])

                stem_vec = (mids[1] - mids[0]).get_array()

                self.assertTrue(allclose(0., dot(stem_vec, twists[0])))
                self.assertTrue(allclose(0., dot(stem_vec, twists[1])))


class TestVectorFunctions(unittest.TestCase):
    '''
    Tests for some of the vector functions in corgy.utilities.vector.
    '''

    def test_change_basis1(self):
        vec1 = array([1,1])

        old_basis = get_standard_basis(len(vec1))
        new_basis = array([[1,1],[1,-1]])

        new_vec = change_basis(vec1, new_basis, old_basis)

        self.assertEqual(new_vec[0], 1)
        self.assertEqual(new_vec[1], 0)

        self.assertTrue(allclose(new_vec, array([1,0])))

    def test_change_basis2(self):
        '''
        Examples from:

        http://tutorial.math.lamar.edu/Classes/LinAlg/ChangeOfBasis.aspx#VS_ChangeBasis_Ex5c
        '''

        vec1 = array([-2., 3., 4.])
        #vec1 = array([1., 1., 1.])

        B = get_standard_basis(len(vec1))
        C = array([[1,-1,1],[0,1,2],[3,0,-1]])

        new_vec = change_basis(vec1, B, C)
        self.assertTrue(allclose(new_vec, array([10,5,0])))

        new_vec = change_basis(array([9,-1,-8]), B, C)
        self.assertTrue(allclose(new_vec, array([-15,-10,15])))

        new_vec = change_basis(array([10,5,0]), C, B)
        self.assertTrue(allclose(new_vec, array([-2.,3.,4.])))

        new_vec = change_basis(array([-6,7,2]), C, B)
        self.assertTrue(allclose(new_vec, array([-21/5., 14/5., -3/5.])))
        #print >>sys.stderr, "new_vec:", new_vec

    def test_spherical_coordinates(self):
        '''
        Test the parameterization of 3D cartesian coordinates to polar
        coordinates.

        See corgy.utilities.vector.spherical_cartesian_to_polar and
            corgy.utilities.vector.spherical_polar_to_cartesian
        '''

        test_cases = [
                [[1., 0., 0.] , [1., pi/2, 0.]],
                [[0., 1., 0.] , [1., pi/2, pi/2]],
                [[0., 2., 0.] , [2., pi/2, pi/2]],
                [[0., 0., 1.] , [1., 0., 0.]],
                [[1., 1., 1.] , [sqrt(3), asin(sqrt(2/3.)), pi/4]],
                [[1., -1., 1.] , [sqrt(3), asin(sqrt(2/3.)), -pi/4]]]

        for case in test_cases:
            polar = spherical_cartesian_to_polar(case[0])
            cartesian = spherical_polar_to_cartesian(case[1])

            self.assertTrue(allclose(polar, case[1]))
            #print "case[0]:", case[0], "cartesian:", cartesian, "case[1]", case[1]
            self.assertTrue(allclose(case[0], cartesian))

    def test_random_spherical_coordinates(self):
        '''
        Use random values to test the conversion between spherical and cartesian coordinates.
        '''

        for i in range(10):
            cart = [uniform(-1., 1.), uniform(-1., 1.), uniform(-1., 1.)]

            polar = spherical_cartesian_to_polar(cart)
            new_cart = spherical_polar_to_cartesian(polar)
            #print "polar:", polar

            self.assertTrue(allclose(cart, new_cart))

if __name__ == '__main__':
    unittest.main() 

