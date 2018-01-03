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
        cg1 = ftmc.CoarseGrainRNA(dotbracket_str = "(((((.....(((((......)))))..)))))")
        self.sm = SpatialModel(cg1)
        self.sm.sample_stats(self.stat_source_limited)
        self.sm.new_traverse_and_build()
        self.mover = fbmov.Mover(self.stat_source_real)
    def test_move_changes_sm(self):
        coords_old = copy.deepcopy(self.sm.bg.coords)
        log.info(self.mover.move(self.sm))
        self.assertNotEqual(self.sm.bg.coords, coords_old, msg="At least one coordinate should have changed significantly.")
    def test_move_and_reset(self):
        for i in range(10):
            coords_old = copy.deepcopy(self.sm.bg.coords)
            log.info(self.mover.move(self.sm))
            self.mover.revert(self.sm)
            self.assertEqual(self.sm.bg.coords, coords_old)

class TestMLSegmentPairMover(TestMoverBaseClassPublicAPI):
    def setUp(self):
        self.stat_source_real = StatStorage("test/fess/data/real.stats")
        cg = ftmc.CoarseGrainRNA(dotbracket_str = "(((((.(((((......)))))....(((((......))))))))))")
        self.sm = SpatialModel(cg)
        self.sm.sample_stats(self.stat_source_real)
        self.sm.new_traverse_and_build()
        self.cutoff = 20
        self.mover = fbmov.MLSegmentPairMover(self.stat_source_real, self.cutoff)
        self.broken_ml, = [ m for m in self.sm.bg.mloop_iterator() if m not in self.sm.bg.mst ]
        self.mover._get_elem = lambda sm: self.broken_ml

    def test_moves_two_ml_segments(self):
        for i in range(3):
            sm = copy.deepcopy(self.sm)
            stats_old = copy.deepcopy(sm.elem_defs)
            self.mover.move(sm)
            changed = set()
            for elem in sm.elem_defs:
                if elem not in stats_old or sm.elem_defs[elem].pdb_name!=stats_old[elem].pdb_name:
                    changed.add(elem)
            self.assertEqual(len(changed),2, msg = "Changed is {}, should be 2 ml segments".format(changed))
    def test_sum_around_junction(self):
        self.mover.move(self.sm)
        sum_stat = statsum_around_junction(self.sm.bg, "m0")
        assert sum_stat.is_similar_to(ftmstat.IDENTITY_STAT, self.cutoff)
    def test_sampl_sum_around_junction(self):
        self.mover.move(self.sm)
        sum_stat = sampled_statsum_around_junction(self.sm, "m2")
        log.warning("Sum stat for test is %s", sum_stat)
        assert sum_stat.is_similar_to(ftmstat.IDENTITY_STAT, self.cutoff)
    def test_sampled_junction_fragment_deviation(self):
        self.mover.move(self.sm)
        for elem in self.sm.bg.mloop_iterator():
            target_stat, = [ s for s in self.sm.bg.get_stats(elem)
                             if s.ang_type == self.sm.bg.get_angle_type(elem, allow_broken = True) ]
            sampled_pdb_name = self.sm.bg.sampled[elem][0]
            for stat in self.stat_source_real.iterate_stats_for(self.sm.bg, elem):
                if stat.pdb_name == sampled_pdb_name:
                    log.info("Asserting similarity for %s (broken = %s): %s =?= %s", elem, self.broken_ml, stat, target_stat)
                    assert stat.is_similar_to(target_stat, self.cutoff)
                    break
    def plot_sampled_junction_fragment_deviation(self):
        self.mover.move(self.sm)
        cg = self.sm.bg
        import matplotlib.pyplot as plt
        for s in it.chain(cg.stem_iterator(), cg.mloop_iterator()):
            if s == self.broken_ml:
                lab = "broken {}".format(s)
            else:
                lab = s
            plt.plot([cg.coords[s][0][0],cg.coords[s][1][0]], [cg.coords[s][0][1],cg.coords[s][1][1]], "o-", label=lab)

        s1_model = place_new_stem(self.sm.stems["s1"], self.sm.elem_defs["s2"], self.sm.elem_defs[self.broken_ml], (1,0))


        error_stat = ftug.get_error_stat(cg, self.broken_ml, "s1", self.sm.elem_defs[self.broken_ml])
        plt.plot([s1_model.mids[0][0], s1_model.mids[1][0]], [s1_model.mids[0][1], s1_model.mids[1][1]], 's-', label="virtual" )

        sum_stat = sampled_statsum_around_junction(self.sm, "m2")
        s1_model = place_new_stem(self.sm.stems["s2"], self.sm.elem_defs["s2"], sum_stat, (1,0))
        plt.plot([s1_model.mids[0][0], s1_model.mids[1][0]], [s1_model.mids[0][1], s1_model.mids[1][1]], 's-', color="blue", label="virtual_around_circle" )

        plt.legend()
        plt.show()
        log.warning("Broken %s, next %s", self.broken_ml, cg.get_next_ml_segment(self.broken_ml))
        assert sum_stat.is_similar_to(ftmstat.IDENTITY_STAT, self.cutoff)
        log.warning("Asserting %s ==  %s", error_stat, ftmstat.IDENTITY_STAT)
        assert error_stat.is_similar_to(ftmstat.IDENTITY_STAT, self.cutoff)

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
