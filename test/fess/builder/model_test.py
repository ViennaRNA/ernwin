import unittest

import forgi.threedee.model.coarse_grain as ftmc
import fess.builder.models as fbm

class TestModel(unittest.TestCase):
    def setUp(self):
        return

    def test_traverse_and_build(self):
        cg = ftmc.CoarseGrainRNA('test/fess/data/1y26.cg')        
        sm = fbm.SpatialModel(cg)

        sm.traverse_and_build()

