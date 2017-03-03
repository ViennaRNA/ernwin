#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input, #pip install future
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)
__metaclass__=type

import unittest
from collections import Counter
import copy

import numpy as np
import numpy.testing as nptest


import forgi.threedee.model.coarse_grain as ftmc
from fess.builder.models import SpatialModel
from fess.builder.stat_container import StatStorage
import fess.builder.move as fbmov
import forgi.threedee.model.stats as ftmstat

from scipy.stats import chisquare

RAND_REPETITIONS = 200
unittest.TestCase.longMessage = True
class TestMoverBaseClassConstruction(unittest.TestCase):
    longMessage = True
    def setUp(self):
        self.stat_source = StatStorage("test/fess/data/test1.stats")
    def test_mover_init(self):
        mover = fbmov.Mover(self.stat_source)
        self.assertIsInstance(mover, fbmov.Mover)
        self.assertIsInstance(mover.stat_source, StatStorage)

class TestMoverBaseClassPrivateMembers(unittest.TestCase):
    longMessage = True
    def setUp(self):        
        self.stat_source_real = StatStorage("test/fess/data/real.stats")
        self.stat_source_limited = StatStorage("test/fess/data/test1.stats")
        cg = ftmc.CoarseGrainRNA(dotbracket_str = "(((((......)))))", seq="GGCGCAAAAAAGCGCC")
        self.sm = SpatialModel(cg)
        self.sm.sample_stats(self.stat_source_limited)
        self.sm.new_traverse_and_build()
        self.mover_realstats = fbmov.Mover(self.stat_source_real)
        self.mover_limitedstats = fbmov.Mover(self.stat_source_limited)

    def test_get_elem_and_stat_real_stats(self):
        elems = Counter()
        for i in range(RAND_REPETITIONS):
            elem, new_stat = self.mover_realstats._get_elem_and_stat(self.sm)
            elems[elem]+=1
            self.assertIn(elem, ["s0", "h0"], msg = "Only elements from SpatialModel can be picked")
            if elem=="s0":
                self.assertIsInstance(new_stat, ftmstat.StemStat)
            else:
                self.assertIsInstance(new_stat, ftmstat.LoopStat)
        print(list(elems.items()))
        c, p = chisquare(list(elems.values()))
        self.assertGreater(p, 0.05, msg = " We do not sample all elements uniformely with a p-value "
            "below 0.2. We sampled with the following frequencies {}. (You might consider repeating "
            "the test if the frequencies look ok)".format(elems.most_common()))
    def test_get_elem_and_stat_limited_stats(self):
        elems = Counter()
        for i in range(RAND_REPETITIONS):
            elem, new_stat = self.mover_limitedstats._get_elem_and_stat(self.sm)
            elems[elem]+=1
            self.assertIn(elem, ["s0", "h0"], msg = "Only elements from SpatialModel can be picked")
            if elem=="s0":
                self.assertEqual(new_stat.pdb_name, "test:s_0")
            else:
                self.assertEqual(new_stat.pdb_name, "test:h_0")
        print(list(elems.items()))
        c, p = chisquare(list(elems.values()))
        self.assertGreater(p, 0.05, msg = " We do not sample all elements uniformely with a p-value "
            "below 0.2. We sampled with the following frequencies {}. (You might consider repeating "
            "the test if the frequencies look ok)".format(elems.most_common()))

    def test_move(self):
        stat = ftmstat.StemStat("stem new:s_0 5 10.388 2.43294047108")
        oldstat = self.sm.elem_defs["s0"]
        self.mover_realstats._prev_stats = {}
        movestring = self.mover_realstats._move(self.sm, "s0", stat)
        self.assertEqual(self.sm.elem_defs["s0"], stat)
        self.assertEqual(self.mover_realstats._prev_stats["s0"], oldstat)
        self.assertEqual(movestring, "s0:test:s_0->new:s_0;")

class TestMoverBaseClassPublicAPI(unittest.TestCase):
    def setUp(self):        
        self.stat_source_real = StatStorage("test/fess/data/real.stats")
        self.stat_source_limited = StatStorage("test/fess/data/test1.stats")
        cg = ftmc.CoarseGrainRNA(dotbracket_str = "(((((.....(((((......)))))..)))))")
        self.sm = SpatialModel(cg)
        self.sm.sample_stats(self.stat_source_limited)
        self.sm.new_traverse_and_build()
        self.mover = fbmov.Mover(self.stat_source_real)
    def test_move_changes_sm(self):
        coords_old = copy.deepcopy(self.sm.bg.coords)
        self.mover.move(self.sm)
        self.assertNotEqual(self.sm.bg.coords, coords_old, msg="At least one coordinate should have changed significantly.")
    def test_move_and_reset(self):
        coords_old = copy.deepcopy(self.sm.bg.coords)
        self.mover.move(self.sm)
        self.mover.revert(self.sm)
        self.assertEqual(self.sm.bg.coords, coords_old)

class TestNMoverPublicAPI(TestMoverBaseClassPublicAPI):
    def setUp(self):
        super(TestNMoverPublicAPI, self).setUp()
        self.mover = fbmov.NElementMover(self.stat_source_real, 2)
        self.mover4 = fbmov.NElementMover(self.stat_source_real, 4)

    def test_move_changes_multiple2(self):
        stats_old = copy.deepcopy(self.sm.elem_defs)
        self.mover.move(self.sm)
        changed = set()
        for elem in self.sm.elem_defs:
            if self.sm.elem_defs[elem].pdb_name!=stats_old[elem].pdb_name:
                changed.add(elem)
        self.assertEqual(len(changed), 2, "Exactely 2 (different) elements should have been changed.")
        
    def test_move_changes_multiple4(self):
        stats_old = copy.deepcopy(self.sm.elem_defs)
        self.assertEqual(self.mover4._n_moves, 4)
        self.mover4.move(self.sm)
        changed = set()
        for elem in self.sm.elem_defs:
            if self.sm.elem_defs[elem].pdb_name!=stats_old[elem].pdb_name:
                changed.add(elem)
        self.assertEqual(len(changed), 4, "Exactely 4 (different) elements should have been changed.")

class TestExhaustiveExplorerPrivateMembers(TestMoverBaseClassPrivateMembers):
    def setUp(self):
        super(TestExhaustiveExplorerPrivateMembers, self).setUp()
        self.mover_realstats = fbmov.ExhaustiveMover(self.stat_source_real, ["s0"], self.sm )
        self.mover_limitedstats = fbmov.ExhaustiveMover(self.stat_source_limited, ["s0"], self.sm )
    def test_get_elem_and_stat_real_stats(self):
        elems = Counter()
        stats = set()
        for i in range(RAND_REPETITIONS):
            elem, new_stat = self.mover_realstats._get_elem_and_stat(self.sm)
            elems[elem]+=1
            self.assertNotIn(new_stat.pdb_name, stats, "ExhaustiveMover should sample every stat only once!")
            stats.add(new_stat.pdb_name)
            self.assertEqual(elem, "s0")
            self.assertIsInstance(new_stat, ftmstat.StemStat)

    def test_get_elem_and_stat_limited_stats(self):
        elems = Counter()
        with self.assertRaises(StopIteration):
            for i in range(12): #More iterations than stats for s0
                self.mover_limitedstats._get_elem_and_stat(self.sm)
                
class TestExhaustiveExplorerPublicAPI(TestMoverBaseClassPublicAPI):
    def setUp(self):
        super(TestExhaustiveExplorerPublicAPI, self).setUp()
        self.mover_realstats = fbmov.ExhaustiveMover(self.stat_source_real, self.sm, ["s0"])
        self.mover_limitedstats = fbmov.ExhaustiveMover(self.stat_source_limited, self.sm, ["s0"])

class TestConnectedElementMoverPublicAPI(TestNMoverPublicAPI):
    def setUp(self):
        super(TestConnectedElementMoverPublicAPI, self).setUp()
        self.mover = fbmov.ConnectedElementMover(self.stat_source_real, 2)
        self.mover4 = fbmov.ConnectedElementMover(self.stat_source_real, 4)

    def test_move_changes_connected(self):            
        for i in range(RAND_REPETITIONS):
            sm = copy.deepcopy(self.sm)
            stats_old = copy.deepcopy(sm.elem_defs)
            self.mover.move(sm)
            changed = set()
            for elem in sm.elem_defs:
                if sm.elem_defs[elem].pdb_name!=stats_old[elem].pdb_name:
                    changed.add(elem)
            print(changed)
            elem1, elem2 = changed
            self.assertIn(elem1, sm.bg.edges[elem2])
            self.assertIn(elem2, sm.bg.edges[elem1])

class TestWholeMLMoverPublicAPI(TestMoverBaseClassPublicAPI):
    def setUp(self):
        self.stat_source_real = StatStorage("test/fess/data/real.stats")
        cg = ftmc.CoarseGrainRNA(dotbracket_str = "((((((((((......))))).....(((((......)))))....(((((......))))))))))")
        self.sm = SpatialModel(cg)
        self.sm.sample_stats(self.stat_source_real)
        self.sm.new_traverse_and_build()
        self.mover = fbmov.WholeMLMover(self.stat_source_real)
    def test_whole_ml_mover_moves_whole_ml(self):
        for i in range(RAND_REPETITIONS):
            sm = copy.deepcopy(self.sm)
            stats_old = copy.deepcopy(sm.elem_defs)
            self.mover.move(sm)
            changed = set()
            for elem in sm.elem_defs:
                if sm.elem_defs[elem].pdb_name!=stats_old[elem].pdb_name:
                    changed.add(elem)
            if any(elem[0]=="m" for elem in changed):
                self.assertEqual(len(changed),3, msg = "Changed is {}, should be 3 ml segments".format(changed))
            else:
                self.assertEqual(len(changed),1)


class TestConvenienceFunctions(unittest.TestCase):
    def setUp(self):
        self.stat_source = StatStorage("test/fess/data/test1.stats")
        cg = ftmc.CoarseGrainRNA(dotbracket_str = "(((((......)))))")
        self.sm = SpatialModel(cg)

    def test_mover_from_string(self):
        mover = fbmov.mover_from_string("Mover", self.stat_source)
        self.assertIsInstance(mover, fbmov.Mover)
        mover = fbmov.mover_from_string("ExhaustiveMover[m0,m1]", self.stat_source, self.sm)
        self.assertIsInstance(mover, fbmov.ExhaustiveMover)
        mover = fbmov.mover_from_string("NElementMover[3]", self.stat_source, self.sm)
        self.assertIsInstance(mover, fbmov.NElementMover)
        self.assertEqual(mover._n_moves, 3)
        mover = fbmov.mover_from_string("ConnectedElementMover[4]", self.stat_source)
        self.assertIsInstance(mover, fbmov.ConnectedElementMover)
        self.assertEqual(mover._n_moves, 4)
        mover = fbmov.mover_from_string("WholeMLMover", self.stat_source)
        self.assertIsInstance(mover, fbmov.WholeMLMover)

