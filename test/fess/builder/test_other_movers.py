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
import fess.builder._other_movers as fbmov



from .move_test import (TestMoverBaseClassConstruction,
                        TestMoverBaseClassPrivateMembers,
                        TestMoverBaseClassPublicAPI,
                        RAND_REPETITIONS)

log = logging.getLogger(__name__)


unittest.TestCase.longMessage = True


class TestNMoverPublicAPI(TestMoverBaseClassPublicAPI):
    def setUp(self):
        super(TestNMoverPublicAPI, self).setUp()
        self.mover = fbmov.NElementMover(2, self.stat_source_real)
        self.mover4 = fbmov.NElementMover(4, self.stat_source_real)

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

class TestMSTchangingMoverPublicAPI(TestMoverBaseClassPublicAPI):
    def setUp(self):
        self.stat_source_real = StatStorage("test/fess/data/real.stats")
        self.stat_source_limited = StatStorage("test/fess/data/test1.stats")
        cg1 = ftmc.CoarseGrainRNA(dotbracket_str = "((((......))((......))))")
        self.sm = SpatialModel(cg1)
        self.sm.sample_stats(self.stat_source_real)
        self.sm.new_traverse_and_build()
        self.mover = fbmov.Mover(self.stat_source_real)
        self.mover = fbmov.MSTchangingMover(self.stat_source_real)
    def test_change_and_reset_m(self):
        assert "m1" in self.sm.bg.mst
        initial_coords = copy.deepcopy(self.sm.bg.coords)
        old_get_elem = self.mover._get_elem
        self.mover._get_elem = lambda *args:"m1"
        log.info(self.mover.move(self.sm))
        self.assertNotEqual(self.sm.bg.coords, initial_coords, msg="At least one coordinate should have changed significantly.")
        self.mover.revert(self.sm)
        self.assertEqual(self.sm.bg.coords, initial_coords)
        self.mover._get_elem = old_get_elem
    def test_mst_is_changed_and_resetted(self):
        assert "m1" in self.sm.bg.mst
        initial_mst = copy.deepcopy(self.sm.bg.mst)
        old_get_elem = self.mover._get_elem
        self.mover._get_elem = lambda *args:"m1"
        log.info(self.mover.move(self.sm))
        self.assertNotEqual(self.sm.bg.mst, initial_mst, msg="The MST should have changed.")
        self.mover.revert(self.sm)
        self.assertEqual(self.sm.bg.mst, initial_mst)
        self.mover._get_elem = old_get_elem


class TestConnectedElementMoverPublicAPI(TestNMoverPublicAPI):
    def setUp(self):
        super(TestConnectedElementMoverPublicAPI, self).setUp()
        self.mover = fbmov.ConnectedElementMover(2, self.stat_source_real)
        self.mover4 = fbmov.ConnectedElementMover(4, self.stat_source_real)

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


class TestExhaustiveExplorerPrivateMembers(TestMoverBaseClassPrivateMembers):
    def setUp(self):
        super(TestExhaustiveExplorerPrivateMembers, self).setUp()
        self.mover_realstats = fbmov.ExhaustiveMover("s0", stat_source=self.stat_source_real, sm=self.sm )
        self.mover_limitedstats = fbmov.ExhaustiveMover("s0", stat_source=self.stat_source_limited, sm=self.sm )
    def test_get_elem_and_stat_real_stats(self):
        elems = Counter()
        stats = set()
        for i in range(RAND_REPETITIONS):
            elem, new_stat = self.mover_realstats._get_elem_and_stat(self.sm)
            elems[elem]+=1
            self.assertNotIn(new_stat.pdb_name, stats,
                             "ExhaustiveMover with single element should sample every stat only once!")
            stats.add(new_stat.pdb_name)
            self.assertEqual(elem, "s0")
            self.assertIsInstance(new_stat, ftmstat.StemStat)

    def test_get_elem_and_stat_limited_stats(self):
        elems = Counter()
        with self.assertRaises(StopIteration):
            for i in range(12): #More iterations than stats for s0
                self.mover_limitedstats._get_elem_and_stat(self.sm)

    def test_iter_moves_many(self):
        cg = ftmc.CoarseGrainRNA(dotbracket_str = "(((((.....)))))", seq="GGCGCAAAAAGCGCC")
        self.sm = SpatialModel(cg)

        mov = fbmov.ExhaustiveMover("s0", "h0", stat_source=self.stat_source_limited,sm=self.sm )
        move_list = list(mov._iter_moves(self.sm))
        self.assertEqual([x[0] for x in move_list], ["s0", "h0", "h0"])
        self.assertEqual([x[1].pdb_name for x in move_list],
                         ["test:s_0", "test:h_1", "test:h_2"])


class TestExhaustiveExplorerPublicAPI(TestMoverBaseClassPublicAPI):
    def setUp(self):
        super(TestExhaustiveExplorerPublicAPI, self).setUp()
        self.mover_realstats = fbmov.ExhaustiveMover("s0", stat_source=self.stat_source_real, sm=self.sm)
        self.mover_limitedstats = fbmov.ExhaustiveMover("s0", stat_source=self.stat_source_limited, sm=self.sm)

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
