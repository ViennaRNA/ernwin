from __future__ import print_function
import os.path as op

import unittest, copy

import fess.builder.energy as fbe

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.stats as ftms
import forgi.utilities.debug as fud

import numpy as np
import fess.builder.models as fbm
from deepdiff import DeepDiff #Use fork from git@github.com:Bernhard10/deepdiff.git !!!

def assertModelsEqual(sm1, sm2, significant_digits=None):
    ON_DEMAND_KEYS=["build_order" ] #Can be None or [] in one model, as they are calculated on demand.
    for key in sm1.__dict__:
        try:
            if sm1.__dict__[key] != sm2.__dict__[key]:
                if not sm2.__dict__[key]  and key in ON_DEMAND_KEYS:
                    continue
                if not sm1.__dict__[key] and key in ON_DEMAND_KEYS:
                    continue
                if key=="bg":
                    assertModelsEqual(sm1.bg, sm2.bg, significant_digits=significant_digits)
                    continue
                if DeepDiff(sm1.__dict__[key], sm2.__dict__[key], significant_digits=significant_digits)=={}:
                    continue
                assert False, "{}.{} not equal: DeepDiff is {}".format(type(sm1), key, DeepDiff(sm1.__dict__[key], sm2.__dict__[key], significant_digits=significant_digits))
        except:
            if not np.array_equal(sm1.__dict__[key], sm2.__dict__[key]):
                if not sm2.__dict__[key]  and key in ON_DEMAND_KEYS:
                    continue
                if not sm1.__dict__[key] and key in ON_DEMAND_KEYS:
                    continue
                if key=="bg":
                    assertModelsEqual(sm1.bg, sm2.bg, significant_digits=significant_digits)
                    continue
                if DeepDiff(sm1.__dict__[key], sm2.__dict__[key], significant_digits=significant_digits)=={}:
                    continue
                assert False, "{}.{} not equal: DeepDiff is {}".format(type(sm1), key, DeepDiff(sm1.__dict__[key], sm2.__dict__[key], significant_digits=significant_digits))

class TestModifyingMST(unittest.TestCase):
    """Test, whether changing the minimum spanning tree works."""
    def setUp(self):
        self.sm = fbm.SpatialModel(ftmc.CoarseGrainRNA('test/fess/data/4way.cg'))
        self.other_sm = fbm.SpatialModel(ftmc.CoarseGrainRNA('test/fess/data/1GID_A.cg'))
    def test_assertModelsEqual_works(self):
        with self.assertRaises(AssertionError):
            assertModelsEqual(self.sm, self.other_sm)
        sm_copy = copy.deepcopy(self.sm)
        assertModelsEqual(self.sm, sm_copy) #Copy of sm (before building) should be equal to sm
        self.sm.load_sampled_elems()
        self.sm.traverse_and_build()
        print("COPY", sm_copy.bg.sampled)
        print("BUILT", self.sm.bg.sampled)
        assertModelsEqual(self.sm, sm_copy, significant_digits=10) #We built sm, but not its copy. Floating point inaccurracies accumulate making only 10 digits significant!
        sm_copy.load_sampled_elems()
        sm_copy.traverse_and_build()
        assertModelsEqual(self.sm, sm_copy) #After building 2 copies independently, they are still the same
        sm_copy2 = copy.deepcopy(self.sm)
        assertModelsEqual(self.sm, sm_copy) #Copy of sm (after building) should be equal to sm
        possible_stats=self.sm.conf_stats.sample_stats(self.sm.bg, "m1")
        new_stat = possible_stats[2]
        old_stat = self.sm.elem_defs["m1"]
        self.sm.elem_defs["m1"] = new_stat
        self.sm.traverse_and_build()
        with assertRaises(AssertionError):
            assertModelsEqual(self.sm, sm_copy) #SM was changed, while sm_copy was not.
        with assertRaises(AssertionError):
            assertModelsEqual(self.sm, sm_copy2) #SM was changed, while sm_copy2 was not.
        self.sm.elem_defs["m1"] = old_stat
        self.sm.traverse_and_build()
        assertModelsEqual(self.sm, sm_copy) #SM has been reset to be like sm_copy.

    def test_load_save_sampled_elems(self):
        sm_copy = copy.deepcopy(self.sm)
        self.sm.load_sampled_elems()
        self.sm.save_sampled_elems()
        assertModelsEqual(self.sm, sm_copy, 9)

class TestModel(unittest.TestCase):
    def setUp(self):
        self.cg = ftmc.CoarseGrainRNA('test/fess/data/1GID_A.cg')
        self.real_stats_fn = 'test/fess/data/real.stats'
        self.filtered_stats_fn = 'test/fess/data/filtered_stats_1gid.csv'

        self.conf_stats = ftms.ConformationStats(self.real_stats_fn)
        self.filtered_stats = ftms.FilteredConformationStats(self.real_stats_fn, 
                                                             self.filtered_stats_fn)
        
        self.sm = fbm.SpatialModel(self.cg, conf_stats = self.conf_stats)
        self.sm.sample_stats()

        return

    def test_filtered_traverse_and_build(self):
        fcs = self.filtered_stats
        ftms.set_conformation_stats(fcs)
        sm = fbm.SpatialModel(self.cg, conf_stats=fcs)
        sm.sample_stats()

        sm.traverse_and_build()
        sm.bg.to_file('temp.cg')

    def test_sample_elems(self):
        sm = fbm.SpatialModel(self.cg, conf_stats = self.conf_stats)

        sm.sample_stats()
    
    def test_get_random_stem_stats(self):
        self.sm.get_random_stem_stats('s0')

    def test_traverse_and_build(self):
        return
        sm = fbm.SpatialModel(self.cg, conf_stats=self.conf_stats)
        sm.sample_stats()

        sm.traverse_and_build()

        cg = ftmc.CoarseGrainRNA('test/fess/data/1ymo_pk.cg')
        sm = fbm.SpatialModel(cg, conf_stats=self.conf_stats)
        sm.sample_stats()
        sm.traverse_and_build()
        sm.bg.to_file('temp1.cg')

        #pseudoknot
        cg = ftmc.CoarseGrainRNA(op.expanduser('~/doarse/4LVV_A/temp.cg'))
        cg = ftmc.from_pdb(op.expanduser('~/doarse/4LVV_A/temp.pdb'),
                           remove_pseudoknots=False)
        sm = fbm.SpatialModel(cg, conf_stats=self.conf_stats)
        sm.sample_stats()
        sm.traverse_and_build()
        #sm.bg.to_file('temp1.cg')

    """
    def test_new_traverse_and_build(self):
        import time, sys

        # the following tests take a long time to execute so they
        # are commented out
        
        time1 = time.time()
        for i in range(10):
            cg = ftmc.CoarseGrainRNA('test/fess/data/4P8Z.cg')
            sm = fbm.SpatialModel(cg)
            sm.sample_stats()

            sm.constraint_energy = fbe.CombinedEnergy([fbe.RoughJunctionClosureEnergy(), 
                                                    fbe.CoarseStemClashEnergy(), 
                                                    fbe.StemVirtualResClashEnergy()])
            sm.traverse_and_build()
        time2 = time.time()
        print >>sys.stderr, "traverse_and_build, elapsed_time:", time2 - time1

        time1 = time.time()
        for i in range(10):
            cg = ftmc.CoarseGrainRNA('test/fess/data/4P8Z.cg')
            sm = fbm.SpatialModel(cg)
            sm.sample_stats()

            sm.constraint_energy = fbe.CombinedEnergy([fbe.RoughJunctionClosureEnergy(), 
                                                    fbe.CoarseStemClashEnergy(), 
                                                    fbe.StemVirtualResClashEnergy()])
            sm.new_traverse_and_build()
        time2 = time.time()

        print >>sys.stderr, "new_traverse_and_build, elapsed time:", time2 - time1"""
        
