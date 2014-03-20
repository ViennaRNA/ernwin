import unittest

from borgy.builder.rmsd import centered_rmsd

from borgy.utilities.vector import get_inter_distances, get_random_vector
from borgy.utilities.vector import rotation_matrix

from random import uniform

from numpy import array, dot, pi, allclose

class TestRmsdFunctions(unittest.TestCase):
    def test_rmsd1(self):
        '''
        Test the rmsd function.
        '''
        crds1 = array([array([1., 0., 0.]), array([0., 1., 0]), array([0., 0., 1.])])
        crds2 = array([array([1., 0., 0.]), array([0., 1., 0]), array([0., 0., 1.])])

        r = centered_rmsd(crds1, crds2)

        self.assertTrue(allclose(r, 0.))

        crds2 = array([array([1., 0., 0.]), array([0., 1., 0]), array([0., 0., 2.])])
        r = centered_rmsd(crds1, crds2)

        self.assertFalse(allclose(r, 0.))
    
    def test_rmsd2(self):
        '''
        Test the optimal superposition function.
        '''

        crds1 = array([array([1., 0., 0.]), array([0., 1., 0]), array([0., 0., 1.])])

        for i in range(10):
            rot = rotation_matrix(get_random_vector(), uniform(0, 2 * pi))
            crds2 = dot( rot, crds1 )

            r = centered_rmsd(crds1, crds2)
            self.assertTrue(allclose(r, 0., atol=1e-7))

        #sup = optimal_superposition(crds1, crds2)
        #crds3 = dot(sup, crds2)

    def test_rmsd3(self):
        '''
        Test the optimal superposition function.
        '''

        crds1 = array([array([1., 0., 0.]), array([0., 1., 0]), array([0., 0., 1.])])
        dists1 = get_inter_distances(crds1)


        for i in range(10):
            translation = get_random_vector(10)
            rot = rotation_matrix(get_random_vector(), uniform(0, 2 * pi))

            crds2_t = crds1 + translation
            crds2_r = array([dot(rot, c) for c in crds2_t])

            dists2 = get_inter_distances(crds2_r)


            r = centered_rmsd(crds1, crds2_r)

            self.assertTrue(allclose(r, 0., atol=1e-7))
