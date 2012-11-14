import unittest
import numpy as np

import corgy.utilities.debug as cud
import corgy.utilities.vector as cuv

from corgy.utilities.vector import normalize
from corgy.utilities.vector import get_double_alignment_matrix, get_standard_basis, change_basis
from corgy.utilities.vector import vec_angle, get_non_colinear_unit_vector
from corgy.utilities.vector import spherical_cartesian_to_polar, spherical_polar_to_cartesian
from corgy.utilities.vector import rotation_matrix, get_inter_distances
from corgy.utilities.vector import get_random_vector, get_random_vector_pair, get_alignment_matrix

from numpy import array, pi, dot, allclose, sqrt, cross
from math import asin

from random import uniform

class TestVectorFunctions(unittest.TestCase):
    '''
    Tests for some of the vector functions in corgy.utilities.vector.
    '''
    def test_line_segment_intersection(self):
        p1 = np.array([1., 0., 0.])
        p2 = np.array([5., 0., 0.])

        p3 = np.array([2., 0., 0.])
        p4 = np.array([2., 5., 0.])

        (i1, i2) = cuv.line_segment_distance(p1, p2, p3, p4)

        self.assertTrue(np.allclose(cuv.magnitude(i2 - i1), 0.))

        (i1, i2) = cuv.line_segment_distance(np.array([1,0,3]),
                                             np.array([-1,0,5]),
                                             np.array([2,0,2]),
                                             np.array([-1,0,3]))

        self.assertTrue(allclose(cuv.magnitude(i2 - i1), 1.414214))

    def test_get_random_vector_pair(self):
        for i in range(10):
            angle = uniform(0, pi)
            vp = get_random_vector_pair(angle)

            self.assertTrue(allclose(angle, vec_angle(vp[0], vp[1])))

    def test_single_alignment(self):
        '''
        Test the alignment of one vector onto another.
        '''
        for i in range(10):
            vec1 = get_random_vector()
            vec2 = get_random_vector()

            mat = get_alignment_matrix(vec1, vec2)

            new_vec = normalize(dot(mat, vec2))
            
            self.assertTrue(allclose(normalize(new_vec), normalize(vec1)))


    def test_double_alignment_rotation(self):
        '''
        Check the function for aligning two sets of vectors.
        '''
        for i in range(10):
            angle = uniform(0, pi)
            #angle = pi / 2

            vp1 = get_random_vector_pair(angle)
            vp2 = get_random_vector_pair(angle)

            mat = get_double_alignment_matrix(vp1, vp2)

            nvp = [0,0]

            nvp[0] = dot(mat, vp2[0])
            nvp[1] = dot(mat, vp2[1])

            self.assertTrue(allclose(normalize(vp1[0]), normalize(nvp[0])))
            self.assertTrue(allclose(normalize(vp1[1]), normalize(nvp[1])))

            pass


    def check_non_colinear_unit_vector(self, vec1):
        ncl = get_non_colinear_unit_vector(vec1)

        comp = cross(vec1, ncl)

        self.assertTrue(not allclose(comp, array([0., 0., 0.])))

        self.assertTrue(allclose(dot(comp, vec1), 0.))
        self.assertTrue(allclose(dot(comp, ncl), 0.))

    def test_get_non_colinear_unit_vector(self):
        vec1 = array([0., -1., 0.])
        self.check_non_colinear_unit_vector(vec1)

        for i in range(10):
            vec1 = array([uniform(-1, 1), uniform(-1, 1), uniform(-1, 1)])
            self.check_non_colinear_unit_vector(vec1)

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
            self.assertTrue(allclose(case[0], cartesian))

    def test_random_spherical_coordinates(self):
        '''
        Use random values to test the conversion between spherical and cartesian coordinates.
        '''

        for i in range(10):
            cart = normalize(array([uniform(-1., 1.), uniform(-1., 1.), uniform(-1., 1.)]))

            polar = spherical_cartesian_to_polar(cart)
            new_cart = spherical_polar_to_cartesian(polar)

            self.assertTrue(allclose(cart, new_cart))

    def test_spherical_rotation(self):
        '''
        Test if we can use the generated spherical coordinates to rotate our vector back to the x-axis.
        '''
        for i in range(10):
            cart = normalize(array([uniform(-1., 1.), uniform(-1., 1.), uniform(-1., 1.)]))
            (r, u, v) = spherical_cartesian_to_polar(cart)


            rot_mat1 = rotation_matrix(array([0., 0., 1.]), v)
            rot_mat2 = rotation_matrix(array([0., 1., 0.]), u - pi/2) 

            new_cart = dot(rot_mat1, cart)
            new_cart = dot(rot_mat2, new_cart)

            self.assertTrue(allclose(new_cart, array([1., 0., 0.])))
