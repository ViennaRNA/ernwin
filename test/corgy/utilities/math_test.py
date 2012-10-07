import unittest, os, math

import corgy.utilities.my_math as cum #hehe
import numpy as np

class TestMath(unittest.TestCase):
    def setUp(self):
        pass

    def test_clock_angle(self):
        self.assertTrue(np.allclose(cum.clock_angle(2., 3.), 1.))
        self.assertTrue(np.allclose(cum.clock_angle(math.pi, 0), math.pi))
        self.assertTrue(np.allclose(cum.clock_angle(math.pi / 2., 0), (3 / 2.) *  math.pi))
