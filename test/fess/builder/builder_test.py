from __future__ import print_function
from __future__ import absolute_import

import tempfile
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, pow, range, round,
                      str, super, zip)

import unittest
import argparse
import copy
import numpy as np
import numpy.testing as nptest
import fess.builder.builder as fbb
import fess.builder.energy as fbe
import fess.builder.models as fbm
from fess.builder import stat_container

import forgi.threedee.model.similarity as ftmsim
import forgi.threedee.model.coarse_grain as ftmc
from six.moves import range

RAND_REPETITION = 10

class TestBuilderAccept(unittest.TestCase):
    def setUp(self):
        cg = ftmc.CoarseGrainRNA.from_bg_file('test/fess/data/4P8Z.cg')
        p = argparse.ArgumentParser()
        fbb.update_parser(p)
        args = p.parse_args([])
        tmpdir = tempfile.mkdtemp(prefix="ernwinTest")
        self.build_function = fbb.from_args(args, None, tmpdir)
        self.sm = fbm.SpatialModel(cg)

    def test_building_rmsd(self):
        orig_cg = copy.deepcopy(self.sm.bg)
        self.build_function(self.sm)
        self.assertAlmostEqual(ftmsim.cg_rmsd(self.sm.bg, orig_cg), 0.)

class TestBuilderBaseClass(unittest.TestCase):
    def setUp(self):
        self.cg = ftmc.CoarseGrainRNA.from_bg_file('test/fess/data/1GID_A-structure1.coord')
        self.cg_copy = ftmc.CoarseGrainRNA.from_bg_file('test/fess/data/1GID_A-structure1.coord')
        self.sm = fbm.SpatialModel(self.cg)

        real_stats_fn = 'fess/stats/all_nr3.36.stats'
        self.stat_source = stat_container.StatStorage(real_stats_fn)
        self.sm.constraint_energy = fbe.StemVirtualResClashEnergy()

        e1 = fbe.RoughJunctionClosureEnergy()
        for ml in self.sm.bg.defines:
            if ml[0]=="m":
                self.sm.junction_constraint_energy[ml] = e1
        self.builder = fbb.Builder(self.stat_source)

    def test_building_samples_stats(self):
        for i in range(RAND_REPETITION):
            self.builder.build(self.sm)
            self.assertGreater(ftmsim.cg_rmsd(self.sm.bg, self.cg_copy), 1)

    def test_clashfree_building(self):
        for i in range(RAND_REPETITION):
            self.builder.build(self.sm)
            self.assertEqual(self.sm.junction_constraint_energy["m0"].eval_energy(self.sm.bg), 0)
            self.assertEqual(self.sm.constraint_energy.eval_energy(self.sm.bg), 0)

class TestChangingMSTBuilder(TestBuilderBaseClass):
    def setUp(self):
        super(TestChangingMSTBuilder, self).setUp()
        self.builder = fbb.ChangingMSTBuilder(self.stat_source)
    def test_building_changes_mst(self):
        # First building
        self.builder.build(self.sm)
        mst = copy.copy(self.sm.bg.mst)
        for i in range(1000):
            self.builder.build(self.sm)
            if self.sm.bg.mst != mst:
                return
        assert False, "The MST was never changed"
