#!/usr/bin/python

import unittest
from sys import stderr

from random import random, uniform

from numpy import array, matrix, allclose, cross
from corgy.utilities.vector import change_basis, get_standard_basis, normalize, create_orthonormal_basis
from corgy.utilities.vector import spherical_cartesian_to_polar, spherical_polar_to_cartesian
from corgy.utilities.vector import rotation_matrix
from corgy.utilities.vector import get_non_colinear_unit_vector
from corgy.utilities.vector import get_vector_centroid
from corgy.utilities.vector import vec_angle

from math import pi, sqrt, asin

from Bio.PDB import PDBParser
from numpy import dot
from corgy.graph.graph_pdb import get_mids, get_twists
from corgy.graph.graph_pdb import get_stem_orientation_parameters
from corgy.graph.bulge_graph import BulgeGraph

from corgy.graph.graph_pdb import get_stem_twist_and_bulge_vecs
from corgy.graph.graph_pdb import get_stem_separation_parameters
from corgy.graph.graph_pdb import stem2_orient_from_stem1, twist2_orient_from_stem1
from corgy.graph.graph_pdb import get_twist_angle, twist2_from_twist1

class TestGraphToAngles(unittest.TestCase):
    '''
    Tests for the gathering of angles statistics.
    '''

    def orientation_test(self, stem1, twist1, stem2, twist2):
        '''
        To test the caculation of the twist angles, the second twist is converted
        to the coordinate system defined by the first stem and its twist.

        The the rotations required for stem2 to be placed on the x-axis are reversed
        and the twist is placed on the x-y plane.
        '''

        (r, u, v, t) = get_stem_orientation_parameters(stem1, twist1, stem2, twist2)
        stem1_basis = create_orthonormal_basis(stem1, twist1)

        twist2_new_basis = change_basis(twist2, stem1_basis, get_standard_basis(3))
        stem2_new_basis = change_basis(stem2, stem1_basis, get_standard_basis(3))

        rot_mat1 = rotation_matrix(array([0., 0., 1.]), v)
        rot_mat2 = rotation_matrix(array([0., 1., 0.]), u - pi/2) 

        twist2_new = dot(rot_mat1, twist2_new_basis)
        twist2_new = dot(rot_mat2, twist2_new)

        stem2_new = dot(rot_mat1, stem2_new_basis)
        stem2_new = dot(rot_mat2, stem2_new)

        #print "stem2:", stem2
        #print "stem2_new:", stem2_new

        self.assertTrue(allclose(dot(stem2_new, array([0., 1., 0.])), 0.))
        self.assertTrue(allclose(dot(stem2_new, array([0., 0., 1.])), 0.))

        self.assertTrue(allclose(dot(twist2_new, array([1., 0., 0.])), 0.))

    def test_get_orientation_parameters(self):
        for i in range(1):
                stem1 = [uniform(-1., 1.), uniform(-1., 1.), uniform(-1., 1.)]
                stem2 = [uniform(-1., 1.), uniform(-1., 1.), uniform(-1., 1.)]

                twist1 = cross(stem1, get_non_colinear_unit_vector(stem1))
                twist2 = cross(stem2, get_non_colinear_unit_vector(stem2))

                self.orientation_test(stem1, twist1, stem2, twist2)
    
    def test_helix_orientations_concrete(self):
        bg = BulgeGraph('test/graph/1gid.comp')

        for define in bg.defines.keys():
            if define[0] != 's' and len(bg.edges[define]) == 2:
                connections = list(bg.edges[define])

                (stem1, twist1, stem2, twist2, bulge) = get_stem_twist_and_bulge_vecs(bg, define, connections)

                self.assertTrue(allclose(dot(stem1, twist1), 0.))
                self.assertTrue(allclose(dot(stem2, twist2), 0.))

    def test_get_helix_parameters(self):
        bg = BulgeGraph('test/graph/1gid.comp')

        for d in bg.defines.keys():
            if d[0] == 's':
                twist_angle = get_twist_angle(bg.coords[d], bg.twists[d])

                self.assertTrue(twist_angle <= pi)
                self.assertTrue(twist_angle >= -pi)

                stem1 = bg.coords[d][1] - bg.coords[d][0]
                twist1 = bg.twists[d][0]

                twist2 = twist2_from_twist1(stem1, twist1, twist_angle)

                self.assertTrue(allclose(dot(stem1, twist2), 0.))
                self.assertTrue(allclose(dot(stem1, twist1), 0.))

                if twist_angle < 0.:
                    twist_angle = -twist_angle

                self.assertTrue(allclose(vec_angle(twist1, twist2), twist_angle))

    def test_reconstruct_helix_orientation(self):
        '''
        Test the reconstruction of a molecule based on the angles.
        '''
        bg = BulgeGraph('test/graph/1gid.comp')

        for define in bg.defines.keys():
            if define[0] != 's' and len(bg.edges[define]) == 2:
                connections = list(bg.edges[define])

                (stem1, twist1, stem2, twist2, bulge) = get_stem_twist_and_bulge_vecs(bg, define, connections)
                (r, u, v, t) = get_stem_orientation_parameters(stem1, twist1, stem2, twist2)
                (r1, u1, v1) = get_stem_separation_parameters(stem1, twist1, bulge)

                stem2_new = stem2_orient_from_stem1(stem1, twist1, (r, u, v))
                self.assertTrue(allclose(stem2_new, stem2))

                twist2_new = twist2_orient_from_stem1(stem1, twist1, (u, v, t))
                self.assertTrue(allclose(twist2_new, twist2))


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

                # make sure the twists are perpendicular to the stem
                self.assertTrue(allclose(0., dot(stem_vec, twists[0])))
                self.assertTrue(allclose(0., dot(stem_vec, twists[1])))

    def test_get_centroid(self):
        '''
        Test the function for getting a centroid.
        '''

        vectors = array([[1,0,0],[2,0,0],[3,0,0]])
        centroid = get_vector_centroid(vectors)

        self.assertTrue(allclose(centroid, array([2, 0, 0])))


class TestVectorFunctions(unittest.TestCase):
    '''
    Tests for some of the vector functions in corgy.utilities.vector.
    '''

    def check_non_colinear_unit_vector(self, vec1):
        ncl = get_non_colinear_unit_vector(vec1)

        comp = cross(vec1, ncl)

        self.assertTrue(not allclose(comp, array([0., 0., 0.])))
        #print "vec:", vec1, "ncl:", ncl, "comp:", comp

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
            cart = normalize(array([uniform(-1., 1.), uniform(-1., 1.), uniform(-1., 1.)]))

            polar = spherical_cartesian_to_polar(cart)
            new_cart = spherical_polar_to_cartesian(polar)
            #print "polar:", polar

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


if __name__ == '__main__':
    unittest.main() 

