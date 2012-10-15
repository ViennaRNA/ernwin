import corgy.utilities.debug as cud

from corgy.graph.bulge_graph import BulgeGraph
from corgy.builder.energy import MissingTargetException 
from corgy.builder.config import Configuration

from corgy.builder.energy import LongRangeInteractionCount
from corgy.builder.energy import RandomEnergy
from corgy.builder.energy import SkewNormalInteractionEnergy
from corgy.builder.energy import CombinedEnergy
from corgy.builder.energy import JunctionClosureEnergy

import corgy.builder.energy as cbe
import corgy.utilities.vector as cuv

import collections

from corgy.builder.models import SpatialModel

import sys, shutil, os, copy

import unittest, pickle

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

    def test_save_and_evaluate_model(self):
        iterations = 12
        output_dir = os.path.join(Configuration.test_output_dir, 'energies')
        energy = pickle.load(open(os.path.join(output_dir, '%s/%d/SkewNormalInteractionEnergy/LongRangeInteractionCount/JunctionClosureEnergy/CombinedEnergy.energy' % (self.bg.name, iterations)), 'r'))

        self.bg.calc_bp_distances()
        sm = SpatialModel(self.bg)
        sm.get_sampled_bulges()
        sm.traverse_and_build()

        sm.bg.output(os.path.join(Configuration.test_output_dir, 'saved_bg.coords'))

        energy1 =  energy.eval_energy(sm, True)

        sm = SpatialModel(BulgeGraph(os.path.join(Configuration.test_output_dir, 'saved_bg.coords')))
        sm.bg.calc_bp_distances()
        sm.get_sampled_bulges()
        #sm.traverse_and_build()
        
        energy2 = energy.eval_energy(sm, True)

        print 'energy1:', energy1
        print 'energy2:', energy2

        self.assertTrue(np.allclose(energy1, energy2))

    def test_stem_clash_energy(self):
        sce = cbe.StemClashEnergy()
        self.sm.build_chain = True
        self.sm.sample_native_stems()
        self.sm.create_native_stem_models()
        energy = sce.eval_energy(self.sm)

        print "energy:", energy

        sm1 = copy.deepcopy(self.sm)
        for i in range(10):
            sm1.sample_stems()
            sm1.sample_angles()
            sm1.traverse_and_build()
            print "energy:", sce.eval_energy(sm1)

    def test_stem_virtual_res_clash_energy(self):
        svrce = cbe.StemVirtualResClashEnergy()
        self.sm.build_chain = True
        self.sm.sample_native_stems()
        self.sm.create_native_stem_models()
        energy = svrce.eval_energy(self.sm)

        print "energy:", energy

class TestDistanceEnergy(unittest.TestCase):
    def test_native_vs_sampled(self):
        '''
        Native energy should always be less than the sampled energy.
        '''
        bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "1y26/graph", "temp.comp"))
        sm = SpatialModel(bg)
        sm.sample_native_stems()
        sm.create_native_stem_models()

        constr = sm.bg.get_long_range_constraints()

        de = cbe.DistanceEnergy(constr)
        energy_native = de.eval_energy(sm)

        self.assertLess(energy_native, 0.2)

        sm.sample_stems()
        sm.sample_angles()
        sm.traverse_and_build()

        energy_sampled = de.eval_energy(sm)
        self.assertLess(energy_native, energy_sampled)

class TestHelixOrientationEnergy(unittest.TestCase):
    def test_eval(self):
        bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))
        sm = SpatialModel(bg)
        sm.sample_native_stems()

        energy = cbe.HelixOrientationEnergy()
        cud.pv('energy.eval_energy(sm)')

        sm.sample_stems()
        sm.sample_angles()
        sm.traverse_and_build()

        energy = cbe.HelixOrientationEnergy()
        cud.pv('energy.eval_energy(sm)')
