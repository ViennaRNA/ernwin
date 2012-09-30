import unittest

import corgy.graph.bulge_graph as cgb
import corgy.visual.force_graph as cvf
import corgy.builder.config as cbc

import os

class TestJSONGraph(unittest.TestCase):
    def setUp(self):
        pass

    def test_convert_graph_to_json(self):
        bg = cgb.BulgeGraph(os.path.join(cbc.Configuration.test_input_dir, "1y26/graph", "temp.comp"))
        #bg = cgb.BulgeGraph(os.path.join(cbc.Configuration.test_input_dir, "1gid/graph", "temp.comp"))

        json_str = cvf.convert_graph_to_json(bg)

        self.assertTrue(len(json_str) > 0)

    def test_convert_graph_to_fancy_json(self):
        bg = cgb.BulgeGraph(os.path.join(cbc.Configuration.test_input_dir, "1y26/graph", "temp.comp"))
        #bg = cgb.BulgeGraph(os.path.join(cbc.Configuration.test_input_dir, "1gid/graph", "temp.comp"))

        print
        json_str = cvf.convert_graph_to_fancy_json(bg)
        print "json_str:", json_str

        self.assertTrue(len(json_str) > 0)
