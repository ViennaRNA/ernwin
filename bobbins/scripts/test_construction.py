#!/usr/bin/python


import unittest
import copy

from corgy.builder.models import AngleStatsDict, AngleStat
from corgy.builder.models import SpatialModel

from corgy.builder.stats import StemStatsDict, LoopStatsDict

from corgy.builder.energy import LongRangeDistanceEnergy, RandomEnergy

from corgy.graph.graph_pdb import get_stem_twist_and_bulge_vecs
from corgy.graph.graph_pdb import get_stem_orientation_parameters
from corgy.graph.graph_pdb import get_stem_separation_parameters

from corgy.graph.bulge_graph import BulgeGraph

from corgy.builder.sampling import GibbsBGSampler, StatisticsPlotter, SamplingStatistics


from math import pi

from random import choice, uniform
from numpy import allclose

class TestStastics(unittest.TestCase):
    '''
    Tests for the statistics loader.
    '''
    def test_loop_stats_loading(self):
        stats_file = '../fess/stats/temp.stats'

        loop_stats = LoopStatsDict(stats_file)

        for dim1 in loop_stats.keys():
            for stat in loop_stats[dim1]:
                self.assertTrue(isinstance(stat.bp_length, ( int, long )))
                self.assertTrue(isinstance(stat.phys_length,  ( float )))

    def test_angle_stats_singleton(self):
        stats_file = '../fess/stats/temp.stats'

        angle_stats1 = AngleStatsDict(stats_file)
        angle_stats2 = AngleStatsDict(stats_file)

        #self.assertEqual(id(angle_stats1), id(angle_stats2))

    def test_angle_stats_loading(self):
        stats_file = '../fess/stats/temp.stats'
    
        angle_stats = AngleStatsDict(stats_file)

        for dim1 in angle_stats.keys():
            for dim2 in angle_stats[dim1].keys():
                for stat in angle_stats[dim1][dim2]:
                    self.assertTrue(isinstance(stat.dim1, ( int, long)))
                    self.assertTrue(isinstance(stat.dim2, ( int, long)))

                    self.assertTrue(isinstance(stat.u, float))
                    self.assertTrue(isinstance(stat.v, float))
                    self.assertTrue(isinstance(stat.t, float))

                    # Polar angles should be in the range [0, pi]
                    self.assertTrue(stat.u <= pi)
                    self.assertTrue(stat.u >= 0.)

                    # Azimuth angles should be in the range [-pi, pi]
                    self.assertTrue(stat.v >= -pi)
                    self.assertTrue(stat.v <= pi)

                    self.assertTrue(isinstance(stat.r1, float))

                    self.assertTrue(isinstance(stat.u1, float))

                    # Polar angles should be in the range [0, pi]
                    self.assertTrue(stat.u1 <= pi)
                    self.assertTrue(stat.u1 >= 0.)

                    self.assertTrue(isinstance(stat.v1, float))
                    
                    # Azimuth angles should be in the range [-pi, pi]
                    self.assertTrue(stat.v1 >= -pi)
                    self.assertTrue(stat.v1 <= pi)

    def test_stem_stats_loading(self):
        stats_file = '../fess/stats/temp.stats'

        stem_stats = StemStatsDict(stats_file)

        for bp_length in stem_stats.keys():
            for stat in stem_stats[bp_length]:

                self.assertTrue(isinstance(stat.bp_length, (int, long)))
                self.assertTrue(isinstance(stat.phys_length, float))
                self.assertTrue(isinstance(stat.twist_angle, float))

                self.assertTrue(stat.twist_angle >= -pi)
                self.assertTrue(stat.twist_angle <= pi)

class TestSpatialModel(unittest.TestCase):


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
        bg = BulgeGraph('../fess/test/graph/1gid.comp')

        sm = SpatialModel(bg)
        sm.traverse_and_build()
        sm.bg.output('this.coords')

        self.check_side_integrity(sm.bg)
        self.check_angle_integrity(sm)

    def test_identical_spatial_model_construction(self):
        bg = BulgeGraph('../fess/test/graph/1gid.comp')

        #self.check_angle_composition(bg, angle_stats)

        sm = SpatialModel(bg)
        sm.traverse_and_build()

        bg1 = copy.deepcopy(sm.bg)
        angle_defs = sm.angle_defs
        stem_defs = sm.stem_defs

        sm1 = SpatialModel(bg, angle_defs = angle_defs, stem_defs = stem_defs)
        sm1.traverse_and_build()
        self.compare_models(bg1, sm1.bg)


    def test_long_range_energy_function(self):
        bg = BulgeGraph('../fess/test/graph/1gid.comp')

        lre = LongRangeDistanceEnergy()
        lre.calibrate(bg, steps=10)
        energy = lre.eval_energy(SpatialModel(bg))

        #print "energy:", energy

        self.assertTrue(isinstance(energy, float))
        self.assertTrue(energy < 0.)


class TestGibbsSampler(unittest.TestCase):

    def test_initiation(self):
        bg = BulgeGraph('../fess/test/graph/1gid.comp')
        bg.calc_bp_distances()
        sm = SpatialModel(bg)

        re = RandomEnergy()

        plotter = StatisticsPlotter()
        random_stats = SamplingStatistics(sm, plotter, 'r', silent=True)
        gs = GibbsBGSampler(sm, re, random_stats)
        gs.step()


if __name__ == '__main__':
    unittest.main()
