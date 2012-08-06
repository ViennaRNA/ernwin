import unittest

import corgy.builder.reconstructor as rc

from corgy.builder.config import Configuration
from corgy.graph.bulge_graph import BulgeGraph

import os

class TestReconstructor(unittest.TestCase):
    def test_reconstruct(self):
        bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "graph", "2b3j.comp"))

        #create random translation, rotation and twist matrices

        #twist, rotate and translate the original stem
        #twist, rotate and translate the model to overlap with the new one


