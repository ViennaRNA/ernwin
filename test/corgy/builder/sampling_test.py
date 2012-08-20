import unittest, os

from corgy.graph.bulge_graph import BulgeGraph
from corgy.builder.config import Configuration
from corgy.builder.models import SpatialModel
from corgy.builder.energy import RandomEnergy
from corgy.builder.sampling import StatisticsPlotter, GibbsBGSampler, SamplingStatistics

class TestGibbsSampler(unittest.TestCase):

    def test_initiation(self):
        bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))
        bg.calc_bp_distances()
        sm = SpatialModel(bg)

        re = RandomEnergy()

        #plotter = StatisticsPlotter()
        plotter = None
        random_stats = SamplingStatistics(sm, plotter, 'r', silent=True)
        gs = GibbsBGSampler(sm, re, random_stats)
        gs.step()
