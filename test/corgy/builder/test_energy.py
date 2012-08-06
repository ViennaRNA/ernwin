from corgy.graph.bulge_graph import BulgeGraph
from corgy.builder.energy import MissingTargetException 

from corgy.builder.energy import LongRangeInteractionCount
from corgy.builder.energy import RandomEnergy
from corgy.builder.energy import SkewNormalInteractionEnergy
from corgy.builder.energy import CombinedEnergy

from corgy.builder.models import SpatialModel

import sys
import unittest

class TestCombinedEnergy(unittest.TestCase):
    def setUp(self):
        self.bg = BulgeGraph('test/files/1gid.comp')
        self.sm = SpatialModel(self.bg)

    def test_calibrate(self):
        snie = SkewNormalInteractionEnergy()
        lric = LongRangeInteractionCount()

        ce = CombinedEnergy([snie, lric])

        ce.calibrate(self.sm)

        lric1 = LongRangeInteractionCount()
        lric1.calibrate(self.sm)

        print "uncalibrated lric:", lric.eval_energy(self.sm)
        print "calibrated lric:", lric1.eval_energy(self.sm)
        
class TestRandomEnergy(unittest.TestCase):
    def setUp(self):
        self.bg = BulgeGraph('test/files/1gid.comp')
        self.sm = SpatialModel(self.bg)
        self.energy_function = RandomEnergy()

    def test_eval_energy(self):
        energy1 = self.energy_function.eval_energy(self.sm)
        energy2 = self.energy_function.eval_energy(self.sm)

        self.assertTrue(energy1 != energy2)

    def test_calibrate(self):
        energy1 = self.energy_function.eval_energy(self.sm)

        self.energy_function.calibrate(self.sm)


class TestLongRangeInteractionCount(unittest.TestCase):
    def setUp(self):
        self.bg = BulgeGraph('test/files/1gid.comp')
        self.sm = SpatialModel(self.bg)

    def test_eval_energy(self):
        lric = LongRangeInteractionCount()

        # the energy function should not be initialized yet
        # so it's target_function variable would not be in existence
        self.assertRaises(MissingTargetException, lric.eval_energy, self.sm)

    def test_calibrate(self):
        lric = LongRangeInteractionCount()

        lric.calibrate(self.sm)
        self.assertTrue(lric.target_interactions != None)
        self.assertTrue(lric.bgf != None)

        print lric.eval_energy(self.sm)

        pass

class TestSkewNormalInteractionEnergy(unittest.TestCase):
    def setUp(self):
        self.bg = BulgeGraph('test/files/1gid.comp')
        self.sm = SpatialModel(self.bg)

    def test_calibrate(self):
        snie = SkewNormalInteractionEnergy()
        snie.calibrate(self.sm, iterations=10)

        pass

