#!/usr/bin/python


import unittest

from corgy.builder.models import AngleStatsDict, AngleStat
from corgy.builder.models import StemStatsDict
from corgy.builder.models import SpatialModel

from corgy.builder.energy import LongRangeDistanceEnergy

from corgy.graph.graph_pdb import get_stem_twist_and_bulge_vecs
from corgy.graph.graph_pdb import get_stem_orientation_parameters
from corgy.graph.graph_pdb import get_stem_separation_parameters

from corgy.graph.bulge_graph import BulgeGraph


from math import pi

from random import choice, uniform
from numpy import allclose

class TestStastics(unittest.TestCase):
    '''
    Tests for the statistics loader.
    '''

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

    def test_spatial_model_construction(self):
        angle_stats = AngleStatsDict('../fess/stats/temp.stats')
        stem_stats = StemStatsDict('../fess/stats/temp.stats')
        bg = BulgeGraph('../fess/test/graph/1gid.comp')

        #self.check_angle_composition(bg, angle_stats)

        sm = SpatialModel(bg, angle_stats, stem_stats)

        sm.traverse_and_build()
        #self.check_angle_composition(sm.bg, angle_stats)

        #self.check_side_integrity(sm.bg)

        sm.bg.output('this.coords')

    def test_long_range_energy_function(self):
        bg = BulgeGraph('../fess/test/graph/1gid.comp')

        lre = LongRangeDistanceEnergy()
        lre.calibrate(bg, steps=10)
        energy = lre.eval_energy(bg)

        #print "energy:", energy

        self.assertTrue(isinstance(energy, float))
        self.assertTrue(energy < 0.)


if __name__ == '__main__':
    unittest.main()
