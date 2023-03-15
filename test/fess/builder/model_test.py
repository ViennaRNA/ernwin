from __future__ import print_function
from __future__ import absolute_import
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, pow, range, round,
                      str, super, zip)

import os.path as op
from pprint import pprint, pformat
import unittest, copy, warnings, sys

import pdb

import fess.builder.energy as fbe

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.similarity as ftmsim
import forgi.threedee.model.stats as ftms
import forgi.threedee.utilities.vector as ftuv
import forgi.utilities.debug as fud
import math
import logging
import numpy as np
import numpy.testing as nptest
import fess.builder.models as fbm
from deepdiff import DeepDiff #Use fork from git@github.com:Bernhard10/deepdiff.git !!!

from fess.builder import stat_container

log=logging.getLogger(__name__)

class TestAddLoop(unittest.TestCase):
    def setUp(self):
        self.example_stem_stat=ftms.StemStat("stem exampleStat 3 5.29399999969 1.19302425058 1 3 7 9")
        self.example_hairpin_stat=ftms.LoopStat("loop exampleStat 3 14.8069260882 1.2124527925 1.12478051025 86 88")
    def test_add_loop_for_hairpin(self):
        cg=ftmc.CoarseGrainRNA.from_dotbracket("(((...)))")
        sm=fbm.SpatialModel(cg)
        sm.elem_defs={}
        sm.elem_defs["s0"]=self.example_stem_stat
        sm.elem_defs["h0"]=self.example_hairpin_stat
        sm.stems['s0'] = sm.add_stem('s0', sm.elem_defs['s0'], fbm.StemModel(),
                                      ftms.AngleStat(), (0,1))
        sm.add_loop("h0","s0")
        self.assertAlmostEqual(ftuv.magnitude(sm.bulges["h0"].mids[1] - sm.bulges["h0"].mids[0]), self.example_hairpin_stat.phys_length)

def value_from_diff(key, obj1, obj2):
    """
    Traverse 1 step in the deepdiff putput hierarchy.

    :param key: e.g. "root.bg.vvecs['s3'][1][0]"
    :param obj1, obj2: The two objects referenced in the key's first part.
                       If key starts with "root", obj1 and obj2 should be of the base type used for comparison.
                       If key starts with bg, they should be CoarseGrainRNAs.
    :returns: A tripple: attribute, remainder, results.
              In the case of a key "root.elem_defs['m1'].pdb_name", the return values would be:
              "elem_defs", "elem_defs['m1'].pdb_name", [obj1.elem_defs['m1'], obj2.elem_defs['m1']]
    """
    attributes=key.split(".")[1:] #([0] is "root")
    remainder=".".join(attributes) #For next iteration
    getitem=[]
    if not attributes:
        return None, None, (None,None)
    attrib=attributes[0]
    if "[" in attrib:
        gets = attrib.split("[")[1:]
        for i, get in enumerate(gets):
            gets[i]=get.partition("]")[0]
    else:
        gets=[]
    attrib=attrib.split("[")[0]
    results=[]
    for a in obj1, obj2:
        try:
            a=a.__dict__[attrib]
            for get in gets:
                if get in a:
                    a=a[get]
        except:
            results.append(None)
        else:
            results.append(a)
    return attrib, remainder, results

def assertModelsEqual(sm1, sm2, significant_digits=None, on_demand_keys=[], ignore_keys=[]):
    """
    :param sm1, sm2: Spatial models that should be compared.
    :significant_digits: INT. Ignore changes in float values more digits than `significant_digits`
                         after the decimal point.
    :on_demand_keys: A list of strings. Add this to ON_DEMAND_KEYS
    :on_demand_keys: A list of strings. Add this to IGNORE_KEYS
    :var ON_DEMAND_KEYS: If this attribute evaluates to False for one of the two spatial models,
                         do not report any changes.
    :var IGNORE_KEYS: Do not report any changes in this attribute.
    :raises: AssertionError, if the models are not equal.
    """
    ON_DEMAND_KEYS=["build_order", "mst", "ang_types", "closed_bulges", "bulges", "newly_added_stems", "stems", "_conf_stats" ]+on_demand_keys
    IGNORE_KEYS = ["newly_added_stems"]+ignore_keys
    diff = DeepDiff(sm1, sm2, significant_digits=significant_digits)
    if diff=={}:
        return
    """
    if "set_item_added" in diff and "set_item_removed" in diff:
        while True:
            for item in  diff["set_item_removed"]:
                if item in diff["set_item_added"]:
                    diff["set_item_removed"].remove(item)
                    diff["set_item_added"].remove(item)
                    break
            else:
                break
    if "iterable_item_added" in diff and "iterable_item_removed" in diff:
        while True:
            for item in  diff["iterable_item_added"]:
                if item in diff["iterable_item_removed"]:
                    del diff["iterable_item_removed"][item]
                    del diff["iterable_item_added"][item]
                    break
            else:
                break"""
    for mode in ["type_changes", "dictionary_item_removed", "dictionary_item_added", "set_item_added", "set_item_removed", "iterable_item_removed", "iterable_item_added", "values_changed"]:
        for key in list(diff.get(mode, {})):
            attribute, remainder, results = value_from_diff(key, sm1, sm2)
            while remainder:
                if (attribute in IGNORE_KEYS or
                   (attribute in ON_DEMAND_KEYS and (not results[0] or not results[1]))):
                    try:
                        del diff[mode][key] #A dictionary
                    except TypeError: #A set
                        diff[mode].remove(key)
                    break
                attribute, remainder, results = value_from_diff(remainder, results[0], results[1])
        if mode in diff and not diff[mode]:
            del diff[mode]
    assert diff=={}, "Deep Diff is\n"+pformat(diff)

def assert_connected_graph(bg):
    build_order = bg.traverse_graph()
    all_stems = set(bg.stem_iterator())
    all_ils = set(bg.iloop_iterator())
    all_mls = set(bg.mloop_iterator())
    built_stems=set()
    for (f, c, t) in build_order:
        built_stems.add(f)
        built_stems.add(t)
        if c in all_mls:
            all_mls.remove(c)
        if c in all_ils:
            all_ils.remove(c)
    error=False
    error_msg=[]
    if len(all_ils) > 0:
        error=True
        error_msg.append("Some interior loops are not part of the build_order: {}".format(all_ils))
    if built_stems!=all_stems:
        error=True
        error_msg.append("Some stems are not part of the build_order: {}".format(all_stems-built_stems))
    for remaining_ml in all_mls:
        connections = bg.connections(remaining_ml)
        if connections[0] not in built_stems:
            error = True
            error_msg.append("Missing multiloop {} connected to missing stem {}".format(remaining_ml, connections[0]))
        if connections[1] not in built_stems:
            error = True
            error_msg.append("Missing multiloop {} connected to missing stem {}".format(remaining_ml, connections[1]))
    assert not error, "\n".join(error_msg)


class TestAsserts(unittest.TestCase):
    """Test, whether the costum asserts defined in this file work as intended."""
    def setUp(self):
        self.sm = fbm.SpatialModel(ftmc.CoarseGrainRNA.from_bg_file('test/fess/data/4way.cg'))
        self.other_sm = fbm.SpatialModel(ftmc.CoarseGrainRNA.from_bg_file('test/fess/data/1GID_A.cg'))
        self.stat_source = stat_container.StatStorage('test/fess/data/test1.stats')
    @unittest.skip("Update for deepdiff V3")
    def test_assertModelsEqual_works(self):
        with self.assertRaises(AssertionError):
            assertModelsEqual(self.sm, self.other_sm)
        sm_copy = copy.deepcopy(self.sm)
        assertModelsEqual(self.sm, sm_copy) #Copy of sm (before building) should be equal to sm
        self.sm.load_sampled_elems()
        self.sm.new_traverse_and_build()
        ignore_if_none=["elem_defs"]
        assertModelsEqual(self.sm, sm_copy, 12, ignore_if_none) #We built sm, but not its copy. Floating point inaccurracies accumulate making only 8 digits significant!
        sm_copy.load_sampled_elems()
        sm_copy.new_traverse_and_build()
        assertModelsEqual(self.sm, sm_copy) #After building 2 copies independently, they are still the same
        sm_copy2 = copy.deepcopy(self.sm)
        assertModelsEqual(self.sm, sm_copy) #Copy of sm (after building) should be equal to sm
        possible_stats=self.stat_source.get_possible_stats(self.sm.bg, "m1")
        new_stat = possible_stats[0]
        old_stat = self.sm.elem_defs["m1"]
        self.sm.elem_defs["m1"] = new_stat
        self.sm.new_traverse_and_build()
        with self.assertRaises(AssertionError):
            assertModelsEqual(self.sm, sm_copy) #SM was changed, while sm_copy was not.
        with self.assertRaises(AssertionError):
            assertModelsEqual(self.sm, sm_copy2) #SM was changed, while sm_copy2 was not.
        self.sm.elem_defs["m1"] = old_stat
        self.sm.new_traverse_and_build()
        assertModelsEqual(self.sm, sm_copy) #SM has been reset to be like sm_copy.
    @unittest.skip("Update for deepdiff V3")
    def test_assertModelsEqual_for_load_save_sampled_elems(self):
        sm_copy = copy.deepcopy(self.sm)
        self.sm.load_sampled_elems()
        self.sm.save_sampled_elems()
        ignore_if_none=["elem_defs"] #elem_defs not loaded in copy.
        assertModelsEqual(sm_copy, self.sm, 9, ignore_if_none)
    def test_assert_connected_graph(self):
        self.sm.bg.traverse_graph()
        assert_connected_graph(self.sm.bg)
        self.sm.bg.mst.remove("m1")
        self.sm.bg.build_order=None
        with self.assertRaises(AssertionError) as e:
            assert_connected_graph(self.sm.bg)
        self.assertIn("Missing multiloop m1", str(e.exception))

class TestStatsFromAndToCoords_IL(unittest.TestCase):
    def setUp(self):
        self.sm1=fbm.SpatialModel(ftmc.CoarseGrainRNA.from_bg_file('test/fess/data/il.cg'))
        self.sm2=fbm.SpatialModel(ftmc.CoarseGrainRNA.from_bg_file('test/fess/data/il.cg'))
    def test_extracting_stats_from_sm_before_building(self):
        self.sm1.load_sampled_elems()
        self.sm2.elem_defs={}
        for d in "i0", "s0", "s1", "h0":
            self.sm2.elem_defs[d] = self.sm1.elem_defs[d]
        self.sm2.new_traverse_and_build()
        self.sm1.new_traverse_and_build()
        self.assertAlmostEqual(ftmsim.cg_rmsd(self.sm1.bg, self.sm2.bg), 0)
        #assertModelsEqual(self.sm1, self.sm2, 12)
    def test_extracting_stats_from_sm_after_building(self):
        self.sm1.load_sampled_elems()
        self.sm1.new_traverse_and_build()
        self.sm2.elem_defs={}
        for d in "i0", "s0", "s1", "h0":
            self.sm2.elem_defs[d] = self.sm1.elem_defs[d]
        self.sm2.new_traverse_and_build()
        self.assertAlmostEqual(ftmsim.cg_rmsd(self.sm1.bg, self.sm2.bg), 0)
        #assertModelsEqual(self.sm1, self.sm2, 12)
    def test_save_and_load_does_not_change_coords(self):
        self.sm1.load_sampled_elems()
        self.sm2.load_sampled_elems()
        self.sm2.new_traverse_and_build()
        cg2_str = self.sm2.bg.to_cg_string()
        bg_loaded=ftmc.CoarseGrainRNA.from_bg_file(cg2_str)
        sm2_reloaded = fbm.SpatialModel( bg_loaded )
        sm2_reloaded.load_sampled_elems()
        self.assertAlmostEqual(ftmsim.cg_rmsd(self.sm2.bg, sm2_reloaded.bg), 0, places=6)
        #assertModelsEqual(self.sm2, sm2_reloaded, 12)
        assertModelsEqual(self.sm1, sm2_reloaded, 12)

class TestStatsFromAndToCoords_ML(unittest.TestCase):
    def setUp(self):
        self.sm1=fbm.SpatialModel(ftmc.CoarseGrainRNA.from_bg_file('test/fess/data/4way.cg'))
        self.sm2=fbm.SpatialModel(ftmc.CoarseGrainRNA.from_bg_file('test/fess/data/4way.cg'))
    def test_extracting_stats_from_sm_before_building(self):
        self.sm1.load_sampled_elems()
        self.sm2.elem_defs={}
        for d in self.sm1.elem_defs:
            self.sm2.elem_defs[d] = self.sm1.elem_defs[d]
        self.sm2.new_traverse_and_build()
        self.sm1.new_traverse_and_build()
        self.assertAlmostEqual(ftmsim.cg_rmsd(self.sm1.bg, self.sm2.bg), 0)
        #assertModelsEqual(self.sm1, self.sm2, 12)
    def test_extracting_stats_from_sm_after_building(self):
        self.sm1.load_sampled_elems()
        self.sm1.new_traverse_and_build()
        self.sm2.elem_defs={}
        for d in self.sm1.elem_defs:
            self.sm2.elem_defs[d] = self.sm1.elem_defs[d]
        self.sm2.new_traverse_and_build()
        self.assertAlmostEqual(ftmsim.cg_rmsd(self.sm1.bg, self.sm2.bg), 0)
        #assertModelsEqual(self.sm1, self.sm2, 12)
    def test_extracting_stats_from_cg_after_building(self):
        self.sm1.load_sampled_elems()
        self.sm1.new_traverse_and_build()
        self.sm2.elem_defs={}
        cg1 = self.sm1.bg
        for d in cg1.defines:
            stats = cg1.get_stats(d)
            self.sm2.elem_defs[d] = stats[0]
        self.sm2.new_traverse_and_build()
        assertModelsEqual(self.sm1, self.sm2, 12, ignore_keys=["sampled", "elem_defs"])
        self.assertAlmostEqual(ftmsim.cg_rmsd(self.sm1.bg, self.sm2.bg), 0)
    def test_building_does_not_change_structure(self):
        self.sm1.load_sampled_elems()
        self.sm2.load_sampled_elems()
        self.sm1.new_traverse_and_build()
        print("RMSD", ftmsim.cg_rmsd(self.sm1.bg, self.sm2.bg))
        assertModelsEqual(self.sm1, self.sm2, 12) #Note: Coords do change!!!
        self.assertAlmostEqual(ftmsim.cg_rmsd(self.sm1.bg, self.sm2.bg), 0, places=6)

    def test_save_and_load_does_not_change_coords(self):
        self.sm1.load_sampled_elems()
        self.sm2.load_sampled_elems()
        self.sm2.new_traverse_and_build()
        cg2_str = self.sm2.bg.to_cg_string()
        bg_loaded=ftmc.CoarseGrainRNA.from_bg_string(cg2_str)
        sm2_reloaded = fbm.SpatialModel( bg_loaded )
        sm2_reloaded.load_sampled_elems()
        assertModelsEqual(self.sm2, sm2_reloaded, 12)
        assertModelsEqual(self.sm1, sm2_reloaded, 12)
        # if py_qcprot is used, an accurracy of only 6 places is expected.
        self.assertAlmostEqual(ftmsim.cg_rmsd(self.sm2.bg, sm2_reloaded.bg), 0, places=6)


class TestModifyingMST(unittest.TestCase):
    """Test, whether changing the minimum spanning tree works."""
    def setUp(self):
        self.sm = fbm.SpatialModel(ftmc.CoarseGrainRNA.from_bg_file('test/fess/data/4way.cg'))
        self.sm.load_sampled_elems()
        self.sm.new_traverse_and_build()
        self.sm_zero_hairpin = fbm.SpatialModel(ftmc.CoarseGrainRNA.from_bg_file('test/fess/data/4GXY_A.cg'))
        self.sm_zero_hairpin.load_sampled_elems()
        self.sm_zero_hairpin.new_traverse_and_build()
        self.sm_pseudoknot = fbm.SpatialModel(ftmc.CoarseGrainRNA.from_bg_file('test/fess/data/pseudoknot.cg'))
        self.sm_pseudoknot.load_sampled_elems()
        self.sm_pseudoknot.new_traverse_and_build()
        self.example_angle_stat=ftms.AngleStat("angle exampleStat 0 1000 1.69462078307 0.313515399557 0.165804917419 5.08692965666 1.04129866007 0.717061903121 3  CC")
    @unittest.skip("Update for deepdiff V3")
    def test_change_and_reset_mst_ML(self):
        self.change_and_reset_mst(self.sm)
    @unittest.skip("Update for deepdiff V3")
    def test_change_and_reset_mst_ZeroLengthHairpin(self):
        self.change_and_reset_mst(self.sm_zero_hairpin)
    @unittest.skip("Update for deepdiff V3")
    def test_how_NOT_to_do_it_get_stats_from_broken_ml_segment_ML(self):
        with self.assertRaises(AssertionError):
            self.how_NOT_to_do_it_get_stats_from_broken_ml_segment(self.sm)
    @unittest.skip("Update for deepdiff V3")
    def test_get_stats_from_broken_ml_segment_ML(self):
        self.get_stats_from_broken_ml_segment(self.sm)
    @unittest.skip("Update for deepdiff V3")
    def test_get_stats_from_broken_ml_segment_ZeroLengthHairpin(self):
        self.get_stats_from_broken_ml_segment(self.sm_zero_hairpin)
    @unittest.skip("Update for deepdiff V3")
    def test_use_sm_set_multiloop_break_segment_ML(self):
        self.use_sm_set_multiloop_break_segment(self.sm)
    @unittest.skip("Update for deepdiff V3")
    def test_use_sm_set_multiloop_break_segment_ZeroLengthHairpin(self):
        self.use_sm_set_multiloop_break_segment(self.sm_zero_hairpin)
    @unittest.skip("Update for deepdiff V3")
    def test_use_sm_set_multiloop_break_segment_Pseudoknot(self):
        self.use_sm_set_multiloop_break_segment(self.sm_pseudoknot)

    def change_and_reset_mst(self, sm):
        sm_copy=copy.deepcopy(sm)
        print("Original Bulges", np.array(sm.bulges["h1"]),"\n",
              ".........Coords", sm.bg.coords["h1"])
        old_stats=sm.elem_defs["m1"]
        junction_nodes = set( x for x in sm.bg.find_bulge_loop("m1", 200) if x[0]=="m" )
        missing_nodes = junction_nodes - sm.bg.mst
        missing_node, = missing_nodes
        #Changing the minimal spanning tree:
        sm.bg.mst.remove("m1")
        sm.bg.mst |= missing_nodes
        sm.bg.build_order = None #No longer valid
        sm.bg.ang_types = None
        del sm.elem_defs["m1"]
        sm.load_sampled_elems()
        print("Load sampled Elems Bulges", np.array(sm.bulges["h1"]),"\n",
              "...................Coords", sm.bg.coords["h1"])
        sm.elem_defs[missing_node] = self.example_angle_stat
        sm.new_traverse_and_build(start='start')
        print("Traversed Bulges", np.array(sm.bulges["h1"]),"\n",
              "..........Coords", sm.bg.coords["h1"])
        #Now change it back to the original state
        sm.bg.mst.remove(missing_node)
        sm.bg.mst |= set(["m1"])
        sm.bg.build_order = None #No longer valid
        sm.bg.ang_types = None
        del sm.elem_defs[missing_node]
        del sm.bg.sampled[missing_node] #We need to delete this!
        sm.load_sampled_elems()
        print("Re-Loaded Original Bulges", np.array(sm.bulges["h1"]),"\n",
              "...................Coords", sm.bg.coords["h1"])
        sm.elem_defs["m1"] = old_stats
        sm.new_traverse_and_build(start='start')
        print("Built again Bulges", np.array(sm.bulges["h1"]),"\n",
              "............Coords", sm.bg.coords["h1"])
        assertModelsEqual(sm_copy, sm, 10)
    def how_NOT_to_do_it_get_stats_from_broken_ml_segment(self, sm):
        """
        This is left in the code to demonstrate, why after changing the minimal spanning tree,
        we must call sm.load_sampled_elems again.
        Assigning the stat to the newly added multiloop segment is NOT enough!!!
        """
        sm_copy=copy.deepcopy(sm)
        junction_nodes = set( x for x in sm.bg.find_bulge_loop("m1", 200) if x[0]=="m" )
        missing_nodes = junction_nodes - sm.bg.mst
        missing_node, = missing_nodes
        #Changing the minimal spanning tree:
        sm.bg.mst.remove("m1")
        sm.bg.mst |= missing_nodes
        sm.bg.build_order = None #No longer valid
        sm.bg.ang_types = None
        build_order = sm.bg.traverse_graph()
        for bo in build_order:
            if bo[1]==missing_node:
                missing_stat=sm.bg.get_bulge_angle_stats_core(missing_node,(bo[0],bo[2]))
                break
        sm.elem_defs[missing_node]=missing_stat
        sm.bg.sampled = dict()
        del sm.elem_defs["m1"]
        sm.new_traverse_and_build(start='start')
        assertModelsEqual(sm_copy, sm, 11, ignore_keys=["build_order", "mst", "sampled", "elem_defs"]) #sm has different mst and the stats for the missing element from the new mst

    def get_stats_from_broken_ml_segment(self, sm):
        sm_copy=copy.deepcopy(sm)
        junction_nodes = set( x for x in sm.bg.find_bulge_loop("m1", 200) if x[0]=="m" )
        missing_nodes = junction_nodes - sm.bg.mst
        missing_node, = missing_nodes
        #Changing the minimal spanning tree:
        sm.bg.mst.remove("m1")
        sm.bg.mst |= missing_nodes
        sm.bg.build_order = None #No longer valid
        sm.bg.ang_types = None
        #del sm.bg.sampled["m1"]
        del sm.elem_defs["m1"]
        build_order = sm.bg.traverse_graph()
        sm.load_sampled_elems()
        sm.new_traverse_and_build(start='start')
        assertModelsEqual(sm_copy, sm, 11, ignore_keys=["build_order", "mst", "sampled", "elem_defs"]) #sm has different mst and the stats for the missing element from the new mst

    def use_sm_set_multiloop_break_segment(self, sm):
        sm_copy=copy.deepcopy(sm)
        assert "m3" in sm_copy.bg.mst #The ML segment we want to break is not yet broken.
        sm.set_multiloop_break_segment("m3")
        assert_connected_graph(sm.bg) #The graph remains connected.
        self.assertNotIn( "m3", sm.bg.mst) #We managed to break the ML at the desired position
        self.assertNotIn( "m3", sm.elem_defs)
        sm.new_traverse_and_build(start="start")
        assertModelsEqual(sm_copy, sm, 11, ignore_keys=["build_order", "mst", "sampled", "elem_defs"]) #sm has different mst and the stats for the missing element from the new mst
class TestModifyingMSTWithoutLoadSampledElements(unittest.TestCase):
    """Additional tests for changing the minimum spanning tree."""
    def setUp(self):
        self.sm = fbm.SpatialModel(ftmc.CoarseGrainRNA.from_bg_file('test/fess/data/4way.cg'))
        self.sm_zero_hairpin = fbm.SpatialModel(ftmc.CoarseGrainRNA.from_bg_file('test/fess/data/4GXY_A.cg'))
        self.sm_pseudoknot = fbm.SpatialModel(ftmc.CoarseGrainRNA.from_bg_file('test/fess/data/pseudoknot.cg'))
        self.example_angle_stat=ftms.AngleStat("angle exampleStat 0 1000 1.69462078307 0.313515399557 0.165804917419 5.08692965666 1.04129866007 0.717061903121 3  CC")
    def test_set_multiloop_break_segment_4way(self):
        self.use_sm_set_multiloop_break_segment_without_loading_sampled_elems(self.sm)
    def test_set_multiloop_break_segment_sm_pseudoknot(self):
        self.use_sm_set_multiloop_break_segment_without_loading_sampled_elems(self.sm_pseudoknot)
    def test_set_multiloop_break_segment_sm_zero_hairpin(self):
        self.use_sm_set_multiloop_break_segment_without_loading_sampled_elems(self.sm_zero_hairpin)
    def use_sm_set_multiloop_break_segment_without_loading_sampled_elems(self, sm):
        sm.set_multiloop_break_segment("m3")
        assert_connected_graph(sm.bg) #The graph remains connected.
        self.assertNotIn( "m3", sm.bg.mst) #We managed to break the ML at the desired position
        # Commented out: With newer versions, the elem_defs can aldo contain info about broken segments.
        # self.assertNotIn( "m3", sm.elem_defs)

class TestModel(unittest.TestCase):
    def setUp(self):
        self.cg = ftmc.CoarseGrainRNA.from_bg_file('test/fess/data/1GID_A.cg')
        self.real_stats_fn = 'test/fess/data/real.stats'
        self.filtered_stats_fn = 'test/fess/data/filtered_stats_1gid.csv'

        self.filtered_stats = ftms.FilteredConformationStats(self.real_stats_fn,
                                                             self.filtered_stats_fn)

        self.stat_source = stat_container.StatStorage(self.real_stats_fn)
        self.sm = fbm.SpatialModel(self.cg)
        self.sm.sample_stats(self.stat_source)

        return

    def test_sample_elems_doesnot_crash(self):
        sm = fbm.SpatialModel(self.cg)

        sm.sample_stats(self.stat_source)

    @unittest.expectedFailure
    def test_get_random_stem_stats(self):
        self.sm.get_random_stem_stats('s0')

class TestNewTraverse(unittest.TestCase):
    def setUp(self):
        self.cg = ftmc.CoarseGrainRNA.from_bg_file('test/fess/data/4GXY_A.cg')
        self.cg_copy = ftmc.CoarseGrainRNA.from_bg_file('test/fess/data/4GXY_A.cg')
        self.sm = fbm.SpatialModel(self.cg)
        real_stats_fn = 'test/fess/data/real.stats'
        self.stat_source = stat_container.StatStorage(real_stats_fn)
        self.cg2 = ftmc.CoarseGrainRNA.from_dotbracket("(((.(((.(((.(((...))).))).))).)))")
        self.cg3 = ftmc.CoarseGrainRNA.from_dotbracket(dotbracket_str = "(((((......)))))", seq="GGCGCAAAAAAGCGCC")

        #(((.(((.(((.(((...))).))).))).)))
        #sssisssisssissshhhsssisssisssisss
        #000211112220333000333022211112000

        self.sm2 = fbm.SpatialModel(self.cg2)
        self.sm2.sample_stats(self.stat_source)

        self.sm3 = fbm.SpatialModel(self.cg3)
        self.sm3.sample_stats(self.stat_source)

    def test_new_traverse_and_build_return_value(self):
        #self.cg2.print_debug()
        nodes = self.sm2.new_traverse_and_build()
        self.assertEqual(nodes, ["s0", "i0", "s1", "i1", "s2", "i2", "s3"])

    def test_new_traverse_and_build_with_start_s_return_value(self):
        self.sm2.new_traverse_and_build()
        nodes = self.sm2.new_traverse_and_build(start = "s2")
        self.assertEqual(nodes, [ "i2", "s3"])
    def test_new_traverse_and_build_with_start_i_return_value(self):
        self.sm2.new_traverse_and_build()
        nodes = self.sm2.new_traverse_and_build(start = "i1")
        self.assertEqual(nodes, [ "i1", "s2", "i2", "s3"])

    def test_new_traverse_and_build_with_step_return_value(self):
        self.sm2.new_traverse_and_build()
        nodes = self.sm2.new_traverse_and_build(max_steps = 2)
        self.assertEqual(nodes, [ "s0", "i0", "s1", "i1", "s2"])

    def test_new_traverse_and_build_with_step_and_start_return_value(self):
        self.sm2.new_traverse_and_build()
        nodes = self.sm2.new_traverse_and_build(start = "s1", max_steps = 1)
        self.assertEqual(nodes, [ "i1", "s2"])

    def test_new_traverse_and_build_no_angles(self):
        nodes = self.sm3.new_traverse_and_build()
        self.assertEqual(nodes, ["s0"])


    def test_new_traverse_and_build(self):
        self.sm.load_sampled_elems()
        self.sm.new_traverse_and_build()
        for k in self.cg_copy.defines.keys():
            log.info("k: %s file: %s, built %s", k, ftuv.magnitude(self.cg_copy.coords.get_direction(k)),
                                                  ftuv.magnitude(self.sm.bg.coords.get_direction(k)))
        self.assertAlmostEqual(ftmsim.cg_rmsd(self.sm.bg, self.cg_copy), 0, places=6)

    def test_new_traverse_and_build_raises_with_start_and_unbuilt_structure(self):
        self.sm.load_sampled_elems()
        with self.assertRaises(ValueError):
            self.sm.new_traverse_and_build(start="s1")

    def test_new_traverse_and_build_start_doesnt_build_before_start(self):
        self.sm.load_sampled_elems()
        #We need to traverse_and_build at least once from the start!
        self.sm.new_traverse_and_build()
        #Change structure at s0
        self.sm.elem_defs["i8"] = self.stat_source.sample_for(self.cg, "i8")

        #Build only part that did not change
        self.sm.new_traverse_and_build(start="s1")
        self.assertAlmostEqual(ftmsim.cg_rmsd(self.sm.bg, self.cg_copy), 0, places=6)
        #Now build everything including the changed s0
        self.sm.new_traverse_and_build()
        self.assertGreater(ftmsim.cg_rmsd(self.sm.bg, self.cg_copy), 9*10**-6)
    def test_new_traverse_and_build_start_really_builds(self):
        self.sm.load_sampled_elems()
        #We need to traverse_and_build at least once from the start!
        self.sm.new_traverse_and_build()
        #Change structure at s6
        self.sm.elem_defs["i6"] = self.stat_source.sample_for(self.cg, "i6")

        #Build part that did change
        self.sm.new_traverse_and_build(start="s2")
        self.assertGreater(ftmsim.cg_rmsd(self.sm.bg, self.cg_copy), 0)

    def test_new_traverse_and_build_steps_doesnt_build_after(self):
        self.sm.load_sampled_elems()
        #We need to traverse_and_build at least once from the start!
        self.sm.new_traverse_and_build()
        self.sm.elem_defs["i6"] = self.stat_source.sample_for(self.cg, "i6")
        self.sm.new_traverse_and_build(max_steps=2)
        self.assertAlmostEqual(ftmsim.cg_rmsd(self.sm.bg, self.cg_copy), 0, places=6)
        self.sm.new_traverse_and_build()
        self.assertGreater(ftmsim.cg_rmsd(self.sm.bg, self.cg_copy), 0)

    def test_new_traverse_and_build_steps_really_builds(self):
        self.sm.load_sampled_elems()
        #We need to traverse_and_build at least once from the start!
        self.sm.new_traverse_and_build()
        print(self.sm.bg.traverse_graph())
        self.sm.elem_defs["i7"] = self.stat_source.sample_for(self.cg, "i7")
        self.sm.new_traverse_and_build(max_steps=5)
        self.assertGreater(ftmsim.cg_rmsd(self.sm.bg, self.cg_copy), 0)

class ReconstructionTests(unittest.TestCase):
    def test_get_stem_rotation_matrix(self):
        stem1 = fbm.StemModel(mids=(np.array([0.,0.,0.]),np.array([0.,0.,10.])), twists=(np.array([0., 1., 0.]),np.array([0., -1., 0.])))
        stem2 = fbm.StemModel(mids=(np.array([0.,0.,0.]),np.array([10.,0.,0.])), twists=(np.array([0., 1., 1.]),np.array([0., -1., -1.])))
        rot_mat = fbm.get_stem_rotation_matrix(stem1, stem2)
        #RotMat to euler-angles: nghiaho.com/?page_id=846
        theta_x = math.atan2(rot_mat[2,1], rot_mat[2,2])
        theta_y = math.atan2(-rot_mat[2,0], math.sqrt(rot_mat[2,1]**2+rot_mat[2,2]**2))
        theta_z = math.atan2(rot_mat[1,0], rot_mat[0,0])
        self.assertTrue(ftuv.is_almost_parallel(np.dot(stem1.vec(), rot_mat), np.array([10,0,0.])))
        self.assertAlmostEqual(theta_x, 0)
        self.assertAlmostEqual(theta_y, -math.pi/2)
        self.assertAlmostEqual(theta_z, -math.pi/4)
