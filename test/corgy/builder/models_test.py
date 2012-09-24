import unittest, os, copy

from corgy.builder.config import Configuration
from corgy.builder.models import SpatialModel
from corgy.builder.stats import AngleStat
from corgy.graph.bulge_graph import BulgeGraph

from numpy import allclose, pi
from random import uniform

import numpy as np
import random, time
import corgy.utilities.vector as cuv
import corgy.builder.models as cbm

class TestSpatialModel(unittest.TestCase):
    def setUp(self):
        self.bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))

    def check_side_integrity(self, bg):
        '''
        Test to make sure that the sides of each stem are
        consistent.
        '''
        for d in bg.defines.keys():
            if d[0] == 's':
                for edge in bg.edges[d]:
                    (s1b, s1e) = bg.get_sides(d, edge)

                    stem_mid = bg.coords[d][s1b]
                    bulge_mids = bg.coords[edge]

                    self.assertTrue(allclose(stem_mid, bulge_mids[0]) or allclose(stem_mid, bulge_mids[1]))

    def check_angle_integrity(self, sm):
        '''
        Test to make sure that the angles match up with the ones
        in the defs.
        '''
        for b in sm.sampled_bulges:
            if len(sm.bg.edges[b]) == 2:
                angle_stat = sm.bg.get_bulge_angle_stats(b)

                self.assertTrue((angle_stat[0] == sm.angle_defs[b]) or (angle_stat[1] == sm.angle_defs[b]))

    def check_angle_composition(self, bg, angle_stats):
        for define in bg.defines.keys():
            if define[0] != 's' and len(bg.edges[define]) == 2:
                connections = list(bg.edges[define])

                (stem1, twist1, stem2, twist2, bulge) = get_stem_twist_and_bulge_vecs(bg, define, connections)

                # Get the orientations for orienting these two stems
                (r, u, v, t) = get_stem_orientation_parameters(stem1, twist1, stem2, twist2)
                (r1, u1, v1) = get_stem_separation_parameters(stem1, twist1, bulge)

                dims = bg.get_bulge_dimensions(define)

                this_stat = AngleStat('', dims[0], dims[1], u, v, t, r1, u1, v1)
                stats_list = angle_stats[dims[0]][dims[1]]

                found = False

                for a_s in stats_list:
                    if a_s == this_stat:

                        found=True
                        break

                self.assertTrue(found)
                         

    def test_angle_stat_equality(self):
        '''
        Check if the equality function for two stats works.
        '''

        for i in range(10):
            u = uniform(0., pi)
            v = uniform(-pi, pi)
            t = uniform(-pi, pi)

            r1 = uniform(0, 20)
            u1 = uniform(0., pi)
            v1 = uniform(-pi, pi)

            a_s1 = AngleStat('', 0, 0, u, v, t, r1, u1, v1)
            a_s2 = AngleStat('', 0, 0, u, v, t, r1, u1, v1)
            a_s3 = AngleStat('', 0, 0, v, u, t, r1, u1, v1)

            self.assertTrue(a_s1 == a_s2)

            self.assertTrue(not (a_s1 == a_s3))


    def compare_models(self, bg1, bg2):
        for d in bg1.defines.keys():
            if d[0] == 's':
                self.assertTrue(allclose(bg1.coords[d], bg2.coords[d]))
                self.assertTrue(allclose(bg1.twists[d][0], bg2.twists[d][0]))
                self.assertTrue(allclose(bg1.twists[d][1], bg2.twists[d][1]))

    def test_spatial_model_construction(self):
        bg = self.bg

        sm = SpatialModel(bg)
        sm.traverse_and_build()
        sm.bg.output('this.coords')

        self.check_side_integrity(sm.bg)
        self.check_angle_integrity(sm)

    def test_identical_spatial_model_construction(self):
        bg = self.bg

        #self.check_angle_composition(bg, angle_stats)

        sm = SpatialModel(bg)
        sm.traverse_and_build()

        bg1 = copy.deepcopy(sm.bg)
        angle_defs = sm.angle_defs
        stem_defs = sm.stem_defs

        sm1 = SpatialModel(bg, angle_defs = angle_defs, stem_defs = stem_defs)
        sm1.traverse_and_build()
        self.compare_models(bg1, sm1.bg)

    def test_sampled_bulges(self):
        bg = BulgeGraph(os.path.join(Configuration.test_input_dir, '1gid/graph/temp.comp'))
        sm = SpatialModel(copy.deepcopy(bg))

        sm.traverse_and_build()
        sb1 = sm.sampled_bulges

        sm.get_sampled_bulges()
        sb2 = sm.sampled_bulges
        #print "sb2:", sb2
        #print "visit_order:", sm.visit_order
        #print "pvsit_order:", sm.prev_visit_order

        self.assertEqual(sb1, sb2)

    def long_bulge_check(self, bg):
        sm = SpatialModel(bg)
        
        sm.traverse_and_build()

        sm.bg.output(os.path.join(Configuration.test_output_dir, 'long_bulges.comp'))
        print 'sm.sampled_bulges:', sm.sampled_bulges

        for key in sm.sampled_bulges:

            if len(bg.edges[key]) != 2:
                continue

            if key not in bg.edges.keys():
                continue
            
            le = list(bg.edges[key])
            e1 = le[0]
            e2 = le[1]

            (s1b, s1e) = sm.bg.get_sides(e1, key)
            (s2b, s2e) = sm.bg.get_sides(e2, key)

            c1 = sm.bg.coords[e1][s1b]
            c2 = sm.bg.coords[e2][s2b]

            self.assertTrue(np.allclose(sm.angle_defs[key].r1, cuv.magnitude(np.array(c1) - np.array(c2))))
            #print 'key:', key, 'dist:', cuv.magnitude(np.array(c1) - np.array(c2))



    def test_long_bulges(self):
        bg1 = BulgeGraph(os.path.join(Configuration.test_input_dir, '1y26/graph/temp.comp'))
        #bg2 = BulgeGraph(os.path.join(Configuration.test_input_dir, '1gid/graph/temp.comp'))

        self.long_bulge_check(bg1)
        #self.long_bulge_check(bg2)

    def are_stem_models_equal(self, sm1, sm2):
        sm1.get_sampled_bulges()

        for i in sm1.visit_order:
            if i[0] != 's':
                continue
            if not sm1.stems[i].__eq__(sm2.stems[i]):
                print "sm1.stems[i]", sm1.stems[i]
                print "sm2.stems[i]", sm2.stems[i]
                return False
        print
        return True

    def test_resample_bulge(self):
        bg = BulgeGraph(os.path.join(Configuration.test_input_dir, '1gid/graph/temp.comp'))
        sm = SpatialModel(copy.deepcopy(bg))

        sm.sample_stems()
        sm.sample_angles()

        sm.traverse_and_build()

        bulge = sm.bg.get_random_bulge()
        dims = sm.bg.get_bulge_dimensions(bulge)
        possible_angles = sm.angle_stats[dims[0]][dims[1]]

        new_angle =  random.choice(possible_angles)
        print "bulge:", bulge

        sm1 = copy.deepcopy(sm)
        sm2 = copy.deepcopy(sm)

        self.assertTrue(self.are_stem_models_equal(sm1, sm2))

        sm1.traverse_and_build()
        sm2.traverse_and_build()

        self.assertTrue(self.are_stem_models_equal(sm1, sm2))

        sm1.angle_defs[bulge] = new_angle
        sm2.angle_defs[bulge] = new_angle

        time1 = time.time()
        for i in range(100):
            sm1.traverse_and_build()
        time2 = time.time()
        print "t1:", time2 - time1

        time1 = time.time()
        for i in range(100):
            sm2.traverse_and_build(start=bulge)
        time2 = time.time()
        print "t2:", time2 - time1

        self.assertTrue(self.are_stem_models_equal(sm1, sm2))

    def check_reconstructed_stems(self, sm, chain, stem_names):
        for stem_name in stem_names:
            stem_def = sm.stem_defs[stem_name]
            bg_stem_def = sm.bg.defines[stem_name]
        
            stem = cbm.define_to_stem_model(chain, bg_stem_def)

            #print "stem:", stem
            #print "sm.stems[stem_name]:", sm.stems[stem_name]

            self.assertEqual(stem, sm.stems[stem_name])

    def test_construct_allatom_stems(self):
        bg = BulgeGraph(os.path.join(Configuration.test_input_dir, '1y26/graph/temp.comp'))
        sm = SpatialModel(bg)

        sm.build_chain = True
        sm.traverse_and_build()

        self.check_reconstructed_stems(sm, sm.chain, sm.stem_defs.keys())
        sm.bg.output('this.coords')

