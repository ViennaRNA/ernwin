import unittest

from corgy.builder.stats import get_loop_stats, get_angle_stats, get_stem_stats

from numpy import pi

class TestStatsFunctions(unittest.TestCase):
    '''
    Tests for the statistics loader.
    '''

    def test_loop_stats_loading(self):
        loop_stats = get_loop_stats()

        for dim1 in loop_stats.keys():
            for stat in loop_stats[dim1]:
                self.assertTrue(isinstance(stat.bp_length, ( int, long )))
                self.assertTrue(isinstance(stat.phys_length,  ( float )))

    def test_angle_stats_singleton(self):
        angle_stats1 = get_angle_stats()
        angle_stats2 = get_angle_stats()

    def test_angle_stats_loading(self):
        angle_stats = get_angle_stats()

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
        stem_stats = get_stem_stats()

        for bp_length in stem_stats.keys():
            for stat in stem_stats[bp_length]:

                self.assertTrue(isinstance(stat.bp_length, (int, long)))
                self.assertTrue(isinstance(stat.phys_length, float))
                self.assertTrue(isinstance(stat.twist_angle, float))

                self.assertTrue(stat.twist_angle >= -pi)
                self.assertTrue(stat.twist_angle <= pi)
            
                self.assertTrue(len(stat.define) == 4)

                self.assertTrue(stat.define[0] < stat.define[1])
                self.assertTrue(stat.define[2] < stat.define[3])

