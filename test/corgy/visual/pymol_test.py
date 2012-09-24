import unittest

from corgy.visual.pymol import PymolPrinter
from corgy.builder.models import SpatialModel
from corgy.graph.bulge_graph import BulgeGraph
from corgy.builder.config import Configuration

import os

class TestPymolPrinter(unittest.TestCase):
    def setUp(self):
        self.pp = PymolPrinter()
        self.bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))
        self.pp.coordinates_to_pymol(self.bg)

    def test_pymol_intro_string(self):
        s = self.pp.pymol_intro_string()

        self.assertNotEquals(len(s), 0)
        self.assertGreaterEqual(s.find('cgo'), 0)
        self.assertGreaterEqual(s.find('cmd'), 0)
        self.assertGreaterEqual(s.find('vfont'), 0)

    def test_pymol_outro_string(self):
        s = self.pp.pymol_outro_string()

        self.assertNotEquals(len(s), 0)
        self.assertGreaterEqual(s.find('load_cgo'), 0)

    def test_pymol_segments_string(self):
        s = self.pp.pymol_segments_string()

        self.assertNotEquals(len(s), 0)
        self.assertGreaterEqual(s.find('CYLINDER'), 0)

    def test_pymol_spheres_string(self):
        s = self.pp.pymol_spheres_string()

        self.assertTrue(s != None)

    def test_pymol_text_string(self):
        s = self.pp.pymol_text_string()

        self.assertTrue(s != None)

    def test_pymol_string(self):
        s = self.pp.pymol_string()

        self.assertTrue(s != None)

    def test_dump_pymol_file(self):
        filename = os.path.join(Configuration.test_output_dir, "1gid_bg_pymol")

        self.pp.dump_pymol_file(filename)
        self.assertTrue(os.path.exists(filename + ".pym"))
        self.assertTrue(os.path.exists(filename + ".pml"))

        f = open(filename + ".pym", 'r')
        lines = f.readlines()
        self.assertGreater(len(lines), 0)
