import unittest, os

from corgy.graph.bulge_graph import BulgeGraph

from corgy.graph.graph_pdb import get_mids, get_twists, get_stem_orientation_parameters
from corgy.graph.graph_pdb import get_stem_twist_and_bulge_vecs, stem2_orient_from_stem1
from corgy.graph.graph_pdb import get_stem_separation_parameters, twist2_orient_from_stem1
from corgy.graph.graph_pdb import get_twist_angle, twist2_from_twist1

from corgy.utilities.vector import get_vector_centroid, get_non_colinear_unit_vector
from corgy.utilities.vector import create_orthonormal_basis, change_basis, get_standard_basis
from corgy.utilities.vector import rotation_matrix, vec_angle

from corgy.builder.config import Configuration

from Bio.PDB import PDBParser

from numpy import allclose, dot, array, cross, pi
from random import uniform

class TestGraphToAngles(unittest.TestCase):
    '''
    Tests for the gathering of angles statistics.
    '''
    def setUp(self):
        self.bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))

    def orientation_subtest(self, stem1, twist1, stem2, twist2):
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


        self.assertTrue(allclose(dot(stem2_new, array([0., 1., 0.])), 0.))
        self.assertTrue(allclose(dot(stem2_new, array([0., 0., 1.])), 0.))

        self.assertTrue(allclose(dot(twist2_new, array([1., 0., 0.])), 0.))

    def test_get_orientation_parameters(self):
        for i in range(1):
                stem1 = array([uniform(-1., 1.), uniform(-1., 1.), uniform(-1., 1.)])
                stem2 = array([uniform(-1., 1.), uniform(-1., 1.), uniform(-1., 1.)])

                twist1 = cross(stem1, get_non_colinear_unit_vector(stem1))
                twist2 = cross(stem2, get_non_colinear_unit_vector(stem2))

                self.orientation_subtest(stem1, twist1, stem2, twist2)
    
    def test_helix_orientations_concrete(self):
        bg = self.bg

        for define in bg.defines.keys():
            if define[0] != 's' and len(bg.edges[define]) == 2:
                connections = list(bg.edges[define])

                (stem1, twist1, stem2, twist2, bulge) = get_stem_twist_and_bulge_vecs(bg, define, connections)

                self.assertTrue(allclose(dot(stem1, twist1), 0.))
                self.assertTrue(allclose(dot(stem2, twist2), 0.))

    def test_get_helix_parameters(self):
        bg = self.bg

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
        bg = self.bg

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

        s = PDBParser().get_structure('test', os.path.join(Configuration.test_input_dir, "1gid/prepare", "temp.pdb"))
        chain = list(s.get_chains())[0]

        bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))

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
