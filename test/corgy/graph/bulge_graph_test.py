import unittest, os

from corgy.graph.bulge_graph import BulgeGraph

from corgy.builder.config import Configuration

from numpy import allclose

class TestBulgeGraph(unittest.TestCase):
    def setUp(self):
        self.bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))

    def test_loop_centroid(self):
        '''
        The centroid of a loop should always be the second coordinate.
        '''
        bg = self.bg

        for d in bg.defines.keys():
            if d[0] != 's' and len(bg.edges[d]) == 1:
                connect = list(bg.edges[d])[0]

                (sb, se) = bg.get_sides(connect, d)
                self.assertTrue(allclose(bg.coords[d][0], bg.coords[connect][sb]))

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

    def t1est_bp_distance(self):
        bg = self.bg
        bg.calc_bp_distances()

        for d1 in bg.defines.keys():
            self.assertTrue(bg.bp_distances[d1][d1] == 0)

        self.assertTrue(bg.bp_distances['b5']['b0'] == 0)
        self.assertTrue(bg.bp_distances['b5']['s6'] == 4)
        self.assertTrue(bg.bp_distances['b5']['x5'] == 3)
        
        self.assertTrue(bg.bp_distances['s5']['s4'] == 2)

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

        self.assertEquals(bg.get_flanking_sequence('b15'), 'UCAACAGAUCUUCUGUUGA')

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

        self.assertEquals(bg.get_flanking_sequence('b0',0), 'ACCAAGUCUCA')


    def t1est_get_define_sequences(self):
        bg = self.bg

        def1 = bg.defines['b11']
        self.assertEquals(bg.seq[def1[0]:def1[1]-1], 'GAAU')

        def1 = bg.defines['b18']
        self.assertEquals(bg.seq[def1[0]:def1[1]], 'UGCAGUUC')
        
        def1 = bg.defines['b15']
        self.assertEquals(bg.seq[def1[0]:def1[1]-1], 'UCU')

        def1 = bg.defines['x4']
        self.assertEquals(bg.seq[def1[2]:def1[3]-1], 'A')

        def1 = bg.defines['b0']
        self.assertEquals(bg.seq[def1[0]:def1[1]-1], 'AAG')

    def t1est_get_indeces_into_flanking(self):
        bg = self.bg

        # CCACGCA
        # )).))))
        # 1111111
        # 0000111
        # 6789012
        (a,b,i1,i2) = bg.get_flanking_handles('x4',1)
        self.assertEquals(i1, 1)
        self.assertEquals(i2, 3)

        # ACCAAGUCUCA
        # (((...(((((
        (a,b,i1,i2) = bg.get_flanking_handles('b0')
        self.assertEquals(a, 36)
        self.assertEquals(b, 40)
        self.assertEquals(i1, 2)
        self.assertEquals(i2, 6)
        
        #UCAACAGAUCUUCUGUUGA
        #((((((((...))))))))
        #0123456789012345678
        (a,b,i1,i2) = bg.get_flanking_handles('b15')
        self.assertEquals(a, 133)
        self.assertEquals(b, 137)
        self.assertEquals(i1, 7)
        self.assertEquals(i2, 11)


