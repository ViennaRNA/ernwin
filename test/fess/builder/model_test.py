from __future__ import print_function
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, pow, range, round,
                      str, super, zip)

import os.path as op
from pprint import pprint, pformat
import unittest, copy, warnings, nose, sys

import fess.builder.energy as fbe

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.comparison as ftme
import forgi.threedee.model.stats as ftms
import forgi.threedee.utilities.vector as ftuv
import forgi.utilities.debug as fud

import numpy as np
import fess.builder.models as fbm
from deepdiff import DeepDiff #Use fork from git@github.com:Bernhard10/deepdiff.git !!!


class TestAddLoop(unittest.TestCase):
    def setUp(self):
        self.example_stem_stat=ftms.StemStat("stem exampleStat 3 5.29399999969 1.19302425058 1 3 7 9")
        self.example_hairpin_stat=ftms.LoopStat("loop exampleStat 3 14.8069260882 1.2124527925 1.12478051025 86 88")
    def test_add_loop_for_hairpin(self):
        cg=ftmc.CoarseGrainRNA()
        cg.from_dotbracket("(((...)))")
        sm=fbm.SpatialModel(cg)
        sm.elem_defs={}
        sm.elem_defs["s0"]=self.example_stem_stat
        sm.elem_defs["h0"]=self.example_hairpin_stat        
        sm.stems['s0'] = sm.add_stem('s0', sm.elem_defs['s0'], fbm.StemModel(), 
                                      ftms.AngleStat(), (0,1))
        sm.add_loop("h0","s0")
        self.assertAlmostEqual(ftuv.magnitude(sm.bulges["h0"].mids[1] - sm.bulges["h0"].mids[0]), self.example_hairpin_stat.phys_length)

def assertModelsEqual(sm1, sm2, significant_digits=None, on_demand_keys=[]):
    ON_DEMAND_KEYS=["build_order", "mst", "ang_types", "closed_bulges", "bulges", "newly_added_stems", "stems", "_conf_stats" ]+on_demand_keys
    IGNORE_KEYS = ["newly_added_stems"]
    diff = DeepDiff(sm1, sm2, significant_digits=significant_digits)
    if diff=={}:
        return
    if "type_changes" in diff:
        for key in list(diff["type_changes"].keys()):
            attribute=key.split(".")[-1]
            if attribute in ON_DEMAND_KEYS:
                if diff["type_changes"][key]["oldtype"] == type(None) or  diff["type_changes"][key]["newtype"] == type(None):
                    del diff["type_changes"][key]
        if diff["type_changes"]=={}:
            del diff["type_changes"]
    if "dic_item_removed" in diff:
        for key in list(diff["dic_item_removed"]):
            attribute=key.split(".")[-1]
            attribute=attribute.partition("[")[0]
            if attribute in ON_DEMAND_KEYS:
                if not sm1.__dict__[attribute] or not sm2.__dict__[attribute]: #Does not work for sm.bg.something #TODO
                    diff["dic_item_removed"].remove(key)
        if not diff["dic_item_removed"]:
            del diff["dic_item_removed"]
    if "iterable_item_removed" in diff:
        for key in list(diff["iterable_item_removed"].keys()):
            attribute=key.split(".")[-1]
            attribute=attribute.partition("[")[0]
            if attribute in ON_DEMAND_KEYS:
                if not sm1.__dict__[attribute] or not sm2.__dict__[attribute]: #Does not work for sm.bg.something #TODO
                    del diff["iterable_item_removed"][key]
        if not diff["iterable_item_removed"]:
            del diff["iterable_item_removed"]
    if "values_changed" in diff:
        for key in list(diff["values_changed"].keys()):
            attribute=key.split(".")[-1]
            attribute=attribute.partition("[")[0]
            print (attribute, "->", key)
            if attribute in IGNORE_KEYS:
                del diff["values_changed"][key]
        if not diff["values_changed"]:
            del diff["values_changed"]
    assert diff=={}, "Deep Diff is\n"+pformat(diff)
  
def assertModelsEqual1(sm1, sm2, significant_digits=None, on_demand_keys=[]):
    if hasattr(sm1, "bg"):
        assert ftme.cg_rmsd(sm1.bg, sm2.bg)<10**-14, "{:.20f} not < 10**-14".format(ftme.cg_rmsd(sm1.bg, sm2.bg))
    ON_DEMAND_KEYS=["build_order", "mst", "ang_types", "closed_bulges", "bulges", "newly_added_stems", "stems", "_conf_stats" ]+on_demand_keys #Can be None or [] in one model, as they are calculated on demand.
    for key in sm1.__dict__:
        ####if key=="coords": continue #The coords can be different if the RMSD is low! (The structure can be rotated/ shifted)
        if key=="newly_added_stems": continue
        try:
            if sm1.__dict__[key] != sm2.__dict__[key]:
                if not sm2.__dict__[key]  and key in ON_DEMAND_KEYS:
                    continue
                if not sm1.__dict__[key] and key in ON_DEMAND_KEYS:
                    continue
                if key=="bg":
                    assertModelsEqual(sm1.bg, sm2.bg, significant_digits, on_demand_keys)
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
                    assertModelsEqual(sm1.bg, sm2.bg, significant_digits, on_demand_keys)
                    continue
                if DeepDiff(sm1.__dict__[key], sm2.__dict__[key], significant_digits=significant_digits)=={}:
                    continue
                assert False, "{}.{} not equal: DeepDiff is {}".format(type(sm1), key, DeepDiff(sm1.__dict__[key], sm2.__dict__[key], significant_digits=significant_digits))

class TestStatsFromAndToCoords_IL(unittest.TestCase):
    def setUp(self):
        self.sm1=fbm.SpatialModel(ftmc.CoarseGrainRNA('test/fess/data/il.cg'))
        self.sm2=fbm.SpatialModel(ftmc.CoarseGrainRNA('test/fess/data/il.cg'))
    def test_extracting_stats_from_sm_before_building(self):
        self.sm1.load_sampled_elems()
        self.sm2.elem_defs={}
        for d in "i0", "s0", "s1", "h0":
            self.sm2.elem_defs[d] = self.sm1.elem_defs[d]
        self.sm2.traverse_and_build()
        self.sm1.traverse_and_build()
        self.assertAlmostEqual(ftme.cg_rmsd(self.sm1.bg, self.sm2.bg), 0)
        assertModelsEqual(self.sm1, self.sm2, 12)
    def test_extracting_stats_from_sm_after_building(self):
        self.sm1.load_sampled_elems()        
        self.sm1.traverse_and_build()
        self.sm2.elem_defs={}
        for d in "i0", "s0", "s1", "h0":
            self.sm2.elem_defs[d] = self.sm1.elem_defs[d]
        self.sm2.traverse_and_build()
        self.assertAlmostEqual(ftme.cg_rmsd(self.sm1.bg, self.sm2.bg), 0)
        assertModelsEqual(self.sm1, self.sm2, 12)
    def test_save_and_load_does_not_change_coords(self):
        self.sm1.load_sampled_elems()
        self.sm2.load_sampled_elems()
        self.sm2.traverse_and_build()
        cg2_str = self.sm2.bg.to_cg_string() 
        bg_loaded=ftmc.CoarseGrainRNA()
        bg_loaded.from_cg_string(cg2_str) 
        sm2_reloaded = fbm.SpatialModel( bg_loaded )
        sm2_reloaded.load_sampled_elems()        
        self.assertAlmostEqual(ftme.cg_rmsd(self.sm2.bg, sm2_reloaded.bg), 0)
        assertModelsEqual(self.sm2, sm2_reloaded, 12)
        assertModelsEqual(self.sm1, sm2_reloaded, 12)

class TestStatsFromAndToCoords_ML(unittest.TestCase):
    def setUp(self):
        self.sm1=fbm.SpatialModel(ftmc.CoarseGrainRNA('test/fess/data/4way.cg'))
        self.sm2=fbm.SpatialModel(ftmc.CoarseGrainRNA('test/fess/data/4way.cg'))
    def test_extracting_stats_from_sm_before_building(self):
        self.sm1.load_sampled_elems()
        self.sm2.elem_defs={}
        for d in self.sm1.elem_defs:
            self.sm2.elem_defs[d] = self.sm1.elem_defs[d]
        self.sm2.traverse_and_build()
        self.sm1.traverse_and_build()
        self.assertAlmostEqual(ftme.cg_rmsd(self.sm1.bg, self.sm2.bg), 0)
        assertModelsEqual(self.sm1, self.sm2, 12)
    def test_extracting_stats_from_sm_after_building(self):
        self.sm1.load_sampled_elems()        
        self.sm1.traverse_and_build()
        self.sm2.elem_defs={}
        for d in self.sm1.elem_defs:
            self.sm2.elem_defs[d] = self.sm1.elem_defs[d]
        self.sm2.traverse_and_build()
        self.assertAlmostEqual(ftme.cg_rmsd(self.sm1.bg, self.sm2.bg), 0)
        assertModelsEqual(self.sm1, self.sm2, 12)
    def test_building_does_not_change_structure(self):
        self.sm1.load_sampled_elems()
        self.sm2.load_sampled_elems()
        self.sm1.traverse_and_build()
        print("RMSD", ftme.cg_rmsd(self.sm1.bg, self.sm2.bg))
        self.assertAlmostEqual(ftme.cg_rmsd(self.sm1.bg, self.sm2.bg), 0)
        assertModelsEqual(self.sm1, self.sm2, 12) #Note: Coords do change!!!
    def test_save_and_load_does_not_change_coords(self):
        self.sm1.load_sampled_elems()
        self.sm2.load_sampled_elems()
        self.sm2.traverse_and_build()
        cg2_str = self.sm2.bg.to_cg_string() 
        bg_loaded=ftmc.CoarseGrainRNA()
        bg_loaded.from_cg_string(cg2_str) 
        sm2_reloaded = fbm.SpatialModel( bg_loaded )
        sm2_reloaded.load_sampled_elems()        
        self.assertAlmostEqual(ftme.cg_rmsd(self.sm2.bg, sm2_reloaded.bg), 0)
        assertModelsEqual(self.sm2, sm2_reloaded, 12)
        assertModelsEqual(self.sm1, sm2_reloaded, 12)

class TestAssertModelEqual(unittest.TestCase):
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
        ignore_if_none=["elem_defs"]
        assertModelsEqual(self.sm, sm_copy, 12, ignore_if_none) #We built sm, but not its copy. Floating point inaccurracies accumulate making only 8 digits significant!
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
        with self.assertRaises(AssertionError):
            assertModelsEqual(self.sm, sm_copy) #SM was changed, while sm_copy was not.
        with self.assertRaises(AssertionError):
            assertModelsEqual(self.sm, sm_copy2) #SM was changed, while sm_copy2 was not.
        self.sm.elem_defs["m1"] = old_stat
        self.sm.traverse_and_build()
        assertModelsEqual(self.sm, sm_copy) #SM has been reset to be like sm_copy.

    def test_load_save_sampled_elems(self):
        sm_copy = copy.deepcopy(self.sm)
        self.sm.load_sampled_elems()
        self.sm.save_sampled_elems()
        ignore_if_none=["elem_defs"] #elem_defs not loaded in copy.
        assertModelsEqual(sm_copy, self.sm, 9, ignore_if_none)

class TestModifyingMST(unittest.TestCase):
    """Test, whether changing the minimum spanning tree works."""
    def setUp(self):
        self.sm = fbm.SpatialModel(ftmc.CoarseGrainRNA('test/fess/data/4way.cg'))
        self.sm.load_sampled_elems()
        self.sm.traverse_and_build()
        self.sm_zero_hairpin = fbm.SpatialModel(ftmc.CoarseGrainRNA('test/fess/data/4GXY_A.cg'))
        self.sm_zero_hairpin.load_sampled_elems()
        self.sm_zero_hairpin.traverse_and_build()
        self.example_angle_stat=ftms.AngleStat("angle exampleStat 0 1000 1.69462078307 0.313515399557 0.165804917419 5.08692965666 1.04129866007 0.717061903121 3  CC")
    def test_change_and_reset_mst_ML(self):
        self.change_and_reset_mst(self.sm)

    def test_change_and_reset_mst_ZeroLengthHairpin(self):
        self.change_and_reset_mst(self.sm_zero_hairpin)

    def change_and_reset_mst(self, sm):

        sm_copy=copy.deepcopy(sm)
        old_stats=sm.elem_defs["m1"]
        junction_nodes = set( x for x in sm.bg.find_bulge_loop("m1", 200) if x[0]=="m" )
        missing_nodes = junction_nodes - sm.bg.mst
        missing_node, = missing_nodes
        #Changing the minimal spanning tree:
        sm.bg.mst.remove("m1")        
        sm.bg.mst |= missing_nodes
        sm.bg.build_order = None #No longer valid
        sm.bg.ang_types = None
        sm.bg.sampled = dict()
        del sm.elem_defs["m1"]
        sm.elem_defs[missing_node] = self.example_angle_stat 
        sm.traverse_and_build(start='start')
        #Now change it back to the original state
        sm.elem_defs["m1"] = old_stats
        sm.bg.mst.remove(missing_node)
        sm.bg.mst |= set(["m1"])
        sm.bg.build_order = None #No longer valid
        sm.bg.ang_types = None
        sm.bg.sampled = dict()
        del sm.elem_defs[missing_node]
        sm.traverse_and_build(start='start')
        assertModelsEqual(sm_copy, sm, 12)



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
        
