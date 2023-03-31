#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input, #pip install future
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from six.moves import range
__metaclass__=type

import unittest
from collections import Counter
import copy
import logging
import numpy as np
import numpy.testing as nptest
import itertools as it

from scipy.stats import chisquare
try:
    from unittest import mock #python3
except ImportError:
    import mock

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.model.stats as ftmstat


from fess.builder.models import SpatialModel, place_new_stem
from fess.builder.stat_container import StatStorage
import fess.builder.move as fbmov




log = logging.getLogger(__name__)


RAND_REPETITIONS = 200

unittest.TestCase.longMessage = True


def sampled_statsum_around_junction(sm, first_elem):
    cg = sm.bg
    elem = first_elem
    stat = sm.elem_defs[elem]
    if stat.ang_type in [-2, -3, 4]:
        stat = -stat
    sum_stat = stat
    elem = cg.get_next_ml_segment(elem)
    i=0
    while elem != first_elem:
        i+=1
        stat = sm.elem_defs[elem]
        if stat.ang_type in [-2, -3, 4]:
            stat = -stat
        sum_stat = ftug.sum_of_stats(sum_stat, stat)
        elem = cg.get_next_ml_segment(elem)
        if i>100:
            assert False
    return sum_stat

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
        cg = ftmc.CoarseGrainRNA.from_dotbracket(dotbracket_str = "(((((......)))))", seq="GGCGCAAAAAAGCGCC")
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
        cg1 = ftmc.CoarseGrainRNA.from_dotbracket(dotbracket_str = "(((((.....(((((......)))))..)))))")
        self.sm = SpatialModel(cg1)
        self.sm.sample_stats(self.stat_source_limited)
        self.sm.new_traverse_and_build()
        self.mover = fbmov.Mover(self.stat_source_real)
    def test_move_changes_sm(self):
        """A move should either change coordinates or the broken ML segment"""
        coords_old = copy.deepcopy(self.sm.bg.coords)
        while True:
            info = self.mover.move(self.sm)
            if info.partition(":")[0] in self.sm.bg.mst:
                break
        self.assertNotEqual(self.sm.bg.coords, coords_old, msg="At least one coordinate should have changed significantly.")
    def test_move_and_reset(self):
        for i in range(10):
            coords_old = copy.deepcopy(self.sm.bg.coords)
            log.info(self.mover.move(self.sm))
            self.mover.revert(self.sm)
            self.assertEqual(self.sm.bg.coords, coords_old)

class TestConvenienceFunctions(unittest.TestCase):
    def setUp(self):
        self.stat_source = StatStorage("test/fess/data/test1.stats")
        cg = ftmc.CoarseGrainRNA.from_dotbracket(dotbracket_str = "(((((......)))))")
        self.sm = SpatialModel(cg)

    def test_integration_with_commandlineparsable(self):
        import fess.builder._other_movers as fbomov # Register additional subclasses of Mover
        mover, = fbmov.Mover.from_string("Mover", stat_source=self.stat_source)
        self.assertIsInstance(mover, fbmov.Mover)
        mover, = fbmov.Mover.from_string("ExhaustiveMover[m0,m1]", stat_source=self.stat_source, sm=self.sm)
        self.assertIsInstance(mover, fbomov.ExhaustiveMover)
        self.assertEqual(list(mover._element_names), ["m0", "m1"])
        mover, = fbmov.Mover.from_string("NElementMover[3]", stat_source=self.stat_source)
        with self.assertRaises(TypeError):
            #Missing kwarg stat_source
            fbmov.Mover.from_string("NElementMover[3]")
        self.assertIsInstance(mover, fbomov.NElementMover)
        self.assertEqual(mover._n_moves, 3)
        mover, = fbmov.Mover.from_string("ConnectedElementMover[4]", stat_source=self.stat_source)
        self.assertIsInstance(mover, fbomov.ConnectedElementMover)
        self.assertEqual(mover._n_moves, 4)
        mover, = fbmov.Mover.from_string("WholeMLMover", stat_source=self.stat_source)
        self.assertIsInstance(mover, fbomov.WholeMLMover)
