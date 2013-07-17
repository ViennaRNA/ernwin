import unittest, os
import itertools as it

import numpy as np

import corgy.builder.config as cbc
import corgy.graph.bulge_graph as cgb
import corgy.utilities.debug as cud

from numpy import allclose

import copy, time

class TestBulgeGraph(unittest.TestCase):
    def setUp(self):
        self.bg = cgb.BulgeGraph(os.path.join(cbc.Configuration.test_input_dir, "1gid/graph", "temp.comp"))

    def test_loop_centroid(self):
        '''
        The centroid of a loop should always be the second coordinate.
        '''
        bg = self.bg

        for d in bg.defines.keys():
            if d[0] != 's' and len(bg.edges[d]) == 1:
                connect = list(bg.edges[d])[0]

                (sb, se) = bg.get_sides(connect, d)
                self.assertTrue(np.allclose(bg.coords[d][0], bg.coords[connect][sb]))

    def test_breadth_fist_traversal(self):
        '''
        Test the breadth-first traversal of a graph.
        '''
        bg = self.bg
        path = bg.breadth_first_traversal()

        # the length of the path should be equal to the number of defines 
        #self.assertEqual(len(bg.defines.keys()), len(path))

        # there should be no duplicates in the path
        path_set = set(path)
        self.assertEqual(len(path), len(path_set))

    def test_calc_vres_distances(self):
        bg = cgb.BulgeGraph(os.path.join(cbc.Configuration.test_input_dir, "1y26/graph", "temp.comp"))
        bg.calc_bp_distances()

        print bg.calc_vres_distance('s0', 4, 's0', 5)
        print bg.calc_vres_distance('s1', 0, 's2', 0)
        print bg.calc_vres_distance('s0', 0, 's2', 0)


    def test_bp_distance(self):
        bg = cgb.BulgeGraph(os.path.join(cbc.Configuration.test_input_dir, "1y26/graph", "temp.comp"))
        bg.calc_bp_distances()

        '''
        for d1 in bg.defines.keys():
            self.assertTrue(bg.bp_distances[d1][d1] == 0)

        for s1 in bg.defines.keys():
            for s2 in bg.defines.keys():
                #print s1, s2, bg.closest_sides[s1][s2], bg.closest_sides[s2][s1]
                if s1[0] == 's' and s2[0] == 's':
                    print s1, s2, bg.closest_sides[s1][s2]
        print '-----------'

        for key1 in bg.bp_distances.keys():
            for key2 in bg.bp_distances[key1]:
                print key1, key2, bg.bp_distances[key1][key2]
        '''

        #print bg.bp_distances
        #cud.pv('bg.closest_sides')

        for (k1, k2) in it.combinations(bg.defines.keys(), 2):
            self.assertEquals(bg.bp_distances[k1][k2], bg.bp_distances[k2][k1])
        '''
        self.assertTrue(bg.bp_distances['b5']['b0'] == 0)
        self.assertTrue(bg.bp_distances['b5']['s6'] == 4)
        self.assertTrue(bg.bp_distances['b5']['x5'] == 3)
        
        self.assertTrue(bg.bp_distances['s5']['s4'] == 2)
        '''

    def test_sequence(self):
        bg = self.bg

        self.assertEquals(bg.seq, "GAAUUGCGGGAAAGGGGUCAACAGCCGUUCAGUACCAAGU"
                 "CUCAGGGGAAACUUUGAGAUGGCCUUGCAAAGGGUAUGGU"
                 "AAUAAGCUGACGGACAUGGUCCUAACCACGCAGCCAAGUC"
                 "CUAAGUCAACAGAUCUUCUGUUGAUAUGGAUGCAGUUC")

    def t1est_get_flanking_sequences(self):
        bg = self.bg

        # Test the 5' end
        f1 = bg.get_flanking_region('b11')
        self.assertEquals(f1[0], 1)
        self.assertEquals(f1[1], 8)
        self.assertEquals(bg.get_flanking_sequence('b11'), 'GAAUUGCG')


        # Test the 3' end
        f1 = bg.get_flanking_region('b18')
        self.assertEquals(f1[0], 148)
        self.assertEquals(f1[1], 158)

        self.assertEquals(bg.get_flanking_sequence('b18'), 'GGAUGCAGUUC')

        # Test a loop
        f1 = bg.get_flanking_region('b15')
        self.assertEquals(f1[0], 126)
        self.assertEquals(f1[1], 144)

        self.assertEquals(bg.get_flanking_sequence('b15'),
                          'UCAACAGAUCUUCUGUUGA')

        # Test a bulge
        f1 = bg.get_flanking_region('x4', 0)
        self.assertEquals(f1[0], 5)
        self.assertEquals(f1[1], 10)

        f1 = bg.get_flanking_region('x4', 1)
        self.assertEquals(f1[0], 106)
        self.assertEquals(f1[1], 112)
        self.assertEquals(bg.get_flanking_sequence('x4', 1), 'CCACGCA')

        # Test a junction
        f1 = bg.get_flanking_region('b0', 0)
        self.assertEquals(f1[0], 34)
        self.assertEquals(f1[1], 44)

        self.assertEquals(bg.get_flanking_sequence('b0', 0), 'ACCAAGUCUCA')

    def t1est_get_define_sequences(self):
        bg = self.bg

        def1 = bg.defines['b11']
        self.assertEquals(bg.seq[def1[0]:def1[1] - 1], 'GAAU')

        def1 = bg.defines['b18']
        self.assertEquals(bg.seq[def1[0]:def1[1]], 'UGCAGUUC')

        def1 = bg.defines['b15']
        self.assertEquals(bg.seq[def1[0]:def1[1] - 1], 'UCU')

        def1 = bg.defines['x4']
        self.assertEquals(bg.seq[def1[2]:def1[3] - 1], 'A')

        def1 = bg.defines['b0']
        self.assertEquals(bg.seq[def1[0]:def1[1] - 1], 'AAG')

    def test_copy(self):
        bg = cgb.BulgeGraph(os.path.join(cbc.Configuration.test_input_dir,
                                         "1gid/graph", "temp.comp"))

        bg1 = copy.deepcopy(bg)
        bg2 = bg.copy()

        self.assertFalse(bg1 is bg2)

        print
        iterations = 1000
        t1 = time.time()
        for i in range(iterations):
            bg1 = copy.deepcopy(bg)
        print "t1:", time.time() - t1

        t1 = time.time()
        for i in range(iterations):
            bg2 = bg.copy()

        print "t2:", time.time() - t1

    def t1est_get_indeces_into_flanking(self):
        bg = self.bg

        # CCACGCA
        # )).))))
        # 1111111
        # 0000111
        # 6789012
        (a, b, i1, i2) = bg.get_flanking_handles('x4', 1)
        self.assertEquals(i1, 1)
        self.assertEquals(i2, 3)

        # ACCAAGUCUCA
        # (((...(((((
        (a, b, i1, i2) = bg.get_flanking_handles('b0')
        self.assertEquals(a, 36)
        self.assertEquals(b, 40)
        self.assertEquals(i1, 2)
        self.assertEquals(i2, 6)

        #UCAACAGAUCUUCUGUUGA
        #((((((((...))))))))
        #0123456789012345678
        (a, b, i1, i2) = bg.get_flanking_handles('b15')
        self.assertEquals(a, 133)
        self.assertEquals(b, 137)
        self.assertEquals(i1, 7)
        self.assertEquals(i2, 11)

    def test_get_stem_angle(self):
        bg = cgb.BulgeGraph(os.path.join(cbc.Configuration.test_input_dir, 
                                         "1gid/graph", "temp.comp"))

        stems = [d for d in bg.defines.keys() if d[0] == 's']

        for i in xrange(len(stems)):
            for j in xrange(i+1, len(stems)):
                print stems[i], stems[j], bg.get_stem_angle(stems[i], stems[j])

    def test_get_twists(self):
        bg = cgb.BulgeGraph(os.path.join(cbc.Configuration.test_input_dir, "1gid/graph", "temp.comp"))
    
        for node in bg.defines.keys():
            self.assertNotEqual(bg.get_twists(node), None)

            if len(bg.edges[node]) == 2:
                print bg.get_twists(node)
                self.assertEqual(len(bg.get_twists(node)), 2)

            if len(bg.edges[node]) == 1:
                print bg.get_twists(node)
                self.assertEqual(len(bg.get_twists(node)), 1)

    def test_from_dotbracket(self):
        fn = os.path.join(cbc.Configuration.test_input_dir, "1gid/prepare", "temp.dotplot")
        bg = cgb.BulgeGraph()
        bg.from_dotbracket_file(fn)
        bg.dump()

if __name__ == '__main__':
    import nose
    nose.run()

