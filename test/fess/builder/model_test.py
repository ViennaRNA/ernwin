import unittest

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.stats as ftms
import forgi.utilities.debug as fud

import fess.builder.models as fbm

class TestModel(unittest.TestCase):
    def setUp(self):
        self.cg = ftmc.CoarseGrainRNA('test/fess/data/1y26.cg')        
        self.conf_stats = ftms.ConformationStats('test/fess/data/real.stats')
        self.sm = fbm.SpatialModel(self.cg, conf_stats = self.conf_stats)

        self.sm.sample_stats()

        return

    def test_traverse_and_build(self):
        pass

    def test_sample_elems(self):
        sm = fbm.SpatialModel(self.cg, conf_stats = self.conf_stats)

        sm.sample_stats()
    
    def test_get_random_stem_stats(self):
        self.sm.get_random_stem_stats('s0')

    def test_traverse_and_build(self):
        sm = fbm.SpatialModel(self.cg, conf_stats=self.conf_stats)
        sm.sample_stats()

        sm.traverse_and_build()
