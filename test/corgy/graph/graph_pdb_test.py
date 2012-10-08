import unittest, os

import corgy.graph.graph_pdb as cgg
import corgy.graph.bulge_graph as cgb

from corgy.graph.graph_pdb import get_mids, get_twists, get_stem_orientation_parameters
from corgy.graph.graph_pdb import get_stem_twist_and_bulge_vecs, stem2_orient_from_stem1
from corgy.graph.graph_pdb import get_stem_separation_parameters, twist2_orient_from_stem1
from corgy.graph.graph_pdb import get_twist_angle, twist2_from_twist1

from corgy.utilities.vector import get_vector_centroid, get_non_colinear_unit_vector
from corgy.utilities.vector import create_orthonormal_basis, change_basis, get_standard_basis
from corgy.utilities.vector import rotation_matrix, vec_angle

from corgy.builder.config import Configuration

import corgy.graph.graph_pdb as cggp
import corgy.utilities.vector as cuv
import random, time
import numpy as np

from Bio.PDB import PDBParser

from numpy import allclose, dot, array, cross, pi
from random import uniform

class TestGraphToAngles(unittest.TestCase):
    '''
    Tests for the gathering of angles statistics.
    '''
    def setUp(self):
        self.bg = cgb.BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))

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

        bg = cgb.BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))

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

    def test_new_stem2_pos_from_stem1(self):
        stem1 = cuv.get_random_vector()
        twist1 = cuv.get_random_vector()

        r = random.uniform(1., 3.)
        u = random.uniform(0, np.pi / 2.)
        v = random.uniform(0, 2 * np.pi)

        stem1_basis = cuv.create_orthonormal_basis(stem1, twist1).transpose()

        print
        t1 = time.time()
        for i in range(100):
            sp1 = cggp.stem2_pos_from_stem1(stem1, twist1, (r,u,v))
        print "t1:", time.time() - t1
        t1 = time.time()
        for i in range(100):
            sp2 = cggp.stem2_pos_from_stem1_1(stem1_basis, (r,u,v))
        print "t2:", time.time() - t1

        print
        print "sp1:", sp1
        print "sp2:", sp2
        self.assertTrue(np.allclose(sp1, sp2))

    def test_new_stem2_orient_from_stem1(self):
        stem1 = cuv.get_random_vector()
        twist1 = cuv.get_random_vector()

        r = random.uniform(1., 3.)
        u = random.uniform(0, np.pi / 2.)
        v = random.uniform(0, 2 * np.pi)

        stem1_basis = cuv.create_orthonormal_basis(stem1, twist1).transpose()

        print
        t1 = time.time()
        for i in range(100):
            sp1 = cggp.stem2_orient_from_stem1(stem1, twist1, (r,u,v))
        print "t1:", time.time() - t1
        t1 = time.time()
        for i in range(100):
            sp2 = cggp.stem2_orient_from_stem1_1(stem1_basis, (r,u,v))
        print "t2:", time.time() - t1

        print
        print "sp1:", sp1
        print "sp2:", sp2
        self.assertTrue(np.allclose(sp1, sp2))

    def test_new_twist2_orient_from_stem1(self):
        stem1 = cuv.get_random_vector()
        twist1 = cuv.get_random_vector()

        r = random.uniform(1., 3.)
        u = random.uniform(0, np.pi / 2.)
        v = random.uniform(0, 2 * np.pi)

        stem1_basis = cuv.create_orthonormal_basis(stem1, twist1).transpose()

        print
        t1 = time.time()
        for i in range(100):
            sp1 = cggp.twist2_orient_from_stem1(stem1, twist1, (r,u,v))
        print "t1:", time.time() - t1
        t1 = time.time()
        for i in range(100):
            sp2 = cggp.twist2_orient_from_stem1_1(stem1_basis, (r,u,v))
        print "t2:", time.time() - t1

        print
        print "sp1:", sp1
        print "sp2:", sp2
        self.assertTrue(np.allclose(sp1, sp2))

    def test_virtual_res_3d_pos(self):
        bg = cgb.BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))

        for d in bg.defines.keys():
            if d[0] == 's':
                stem = d
                stem_len = bg.defines[d][1] - bg.defines[d][0] + 1
                (pos, vec) = cgg.virtual_res_3d_pos(bg, stem, 0)

                self.assertTrue(allclose(bg.twists[d][0] + bg.coords[d][0], pos + vec))
                (pos, vec) = cgg.virtual_res_3d_pos(bg, stem, stem_len - 1)

                self.assertTrue(allclose(bg.twists[d][1] + bg.coords[d][1], pos + vec))

    def test_virtual_res_basis(self):
        bg = cgb.BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))

        for d in bg.defines.keys():
            if d[0] == 's':
                stem_len = bg.defines[d][1] - bg.defines[d][0] + 1
                for i in range(stem_len):
                    b1 = cgg.virtual_res_basis(bg, d, i)

    def test_surrounding_pos_to_cartesian(self):
        bg = cgb.BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))
        stems = [d for d in bg.defines.keys() if d[0] == 's']

        for i in range(len(stems)):
            s1_len = bg.defines[stems[i]][1] - bg.defines[stems[i]][0] + 1
            
            for j in range(i+1, len(stems)):
                s2_len = bg.defines[stems[j]][1] - bg.defines[stems[j]][0] + 1

                for k in range(s1_len):
                    for l in range(s2_len):
                        # re-define the position of (stems[j], l) into the 
                        # coordinate system defined by (stems[i], k) 
                        spos = cgg.pos_to_spos(bg, stems[i], k,
                                               stems[j], l)

                        # get the actual position of (stems[j], l) in cartesian space
                        (vpos1, vvec1) = cgg.virtual_res_3d_pos(bg, stems[j], l)
                        vpos1 = vpos1 + vvec1

                        # convert the spos back to the real coordinate system and see
                        # make sure it matches the one calculated above
                        r_vpos1 = spos_to_cartesian(bg, stems[i], k, spos)
                        self.assertTrue(allclose(vpos1, r_vpos1))

