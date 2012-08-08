from corgy.graph.bulge_graph import BulgeGraph
from corgy.builder.energy import MissingTargetException 
from corgy.builder.config import Configuration

from corgy.builder.energy import LongRangeInteractionCount
from corgy.builder.energy import RandomEnergy
from corgy.builder.energy import SkewNormalInteractionEnergy
from corgy.builder.energy import CombinedEnergy
from corgy.builder.energy import JunctionClosureEnergy

from corgy.builder.models import SpatialModel

import sys, shutil, os
import unittest
import numpy as np
        
class TestRandomEnergy(unittest.TestCase):
    def setUp(self):
        self.bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "1y26/graph", "temp.comp"))
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
        self.bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "1y26/graph", "temp.comp"))
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
        self.bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "1y26/graph", "temp.comp"))
        self.sm = SpatialModel(self.bg)

    def test_calibrate(self):
        snie = SkewNormalInteractionEnergy()
        snie.calibrate(self.sm, iterations=10)

        pass

class TestJunctionClosureEnergy(unittest.TestCase):
    def setUp(self):
        self.bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "1y26/graph", "temp.comp"))
        self.sm = SpatialModel(self.bg)

    def test_calibrate(self):
        jce = JunctionClosureEnergy()
        jce.calibrate(self.sm, iterations=10)


class TestCombinedEnergy(unittest.TestCase):
    def setUp(self):
        self.bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "1y26/graph", "temp.comp"))
        self.sm = SpatialModel(self.bg)

    def test_calibrate(self):
        iterations = 12

        snie = SkewNormalInteractionEnergy()
        lric = LongRangeInteractionCount()
        jce = JunctionClosureEnergy()

        output_dir = os.path.join(Configuration.test_output_dir, 'energies')
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        ce = CombinedEnergy([snie, lric, jce])
        ce.calibrate(self.sm, iterations=12, bg_energy = None, output_dir=output_dir)

        lric1 = LongRangeInteractionCount()
        lric1.calibrate(self.sm)

        self.assertTrue(os.path.exists(os.path.join(output_dir, '%s/%d/SkewNormalInteractionEnergy.energy' % (self.bg.name, iterations))))
        self.assertTrue(os.path.exists(os.path.join(output_dir, '%s/%d/SkewNormalInteractionEnergy/LongRangeInteractionCount.energy' % (self.bg.name, iterations))))
        self.assertTrue(os.path.exists(os.path.join(output_dir, '%s/%d/SkewNormalInteractionEnergy/LongRangeInteractionCount/JunctionClosureEnergy.energy' % (self.bg.name, iterations))))
        self.assertTrue(os.path.exists(os.path.join(output_dir, '%s/%d/SkewNormalInteractionEnergy/LongRangeInteractionCount/JunctionClosureEnergy/CombinedEnergy.energy' % (self.bg.name, iterations))))

        energy1 = lric.eval_energy(self.sm)
        energy2 = lric1.eval_energy(self.sm)

        print "uncalibrated lric:", energy1
        print "calibrated lric:", energy2

        self.assertFalse(np.allclose(energy1, energy2, atol=0.1))
