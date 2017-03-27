from __future__ import print_function
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, pow, range, round,
                      str, super, zip)

import unittest
import itertools as it
import numpy as np
import numpy.testing as nptest
import fess.builder.aminor as fba
import forgi.threedee.model.coarse_grain as ftmc
import math
from StringIO import StringIO
try:
    from unittest.mock import mock_open, patch #Python 3
except:
    from mock import mock_open, patch #Python 2

A_MULTICHAIN_FR3D_OUTPUT = StringIO(""" 
Filename Discrepancy       1       2       3 Cha   1-2   1-3   2-3 Con   1-2   1-3   2-3   1-2   1-3   2-1   2-3   3-1   3-2   1-2   1-3   2-1   2-3   3-1   3-2   1-2   1-3   2-3
    1S72      0.0340  A  104  A  957  U 1009 900 ----  ----   cWW  AAA  1909  1885    24                                                                             -     -     -
""")
A_FR3D_OUTPUT = StringIO("""
Filename Discrepancy       1       2       3 Cha   1-2   1-3   2-3 Con   1-2   1-3   2-3   1-2   1-3   2-1   2-3   3-1   3-2   1-2   1-3   2-1   2-3   3-1   3-2   1-2   1-3   2-3
    1S72      0.4119  A  187  C  154  G  182 000 ----  ----   cWW  AAA    33     5    28                                                                             -     -     -
""")

ANOTHER_FR3D_OUTPUT = StringIO("""
Filename Discrepancy       1       2       3 Cha   1-2   1-3   2-3 Con   1-2   1-3   2-3   1-2   1-3   2-1   2-3   3-1   3-2   1-2   1-3   2-1   2-3   3-1   3-2   1-2   1-3   2-3
    1S72      0.4119  A  187  C  154  G  182 000 ----  ----   cWW  AAA    33     5    28                                                                             -     -     -
    1S72      0.4247  A 1476  G 1867  C 1862 000  ---- ----   cWW  AAA   390   385     5                                     n2BR                                    -     -     -
    1S72      0.4416  A  395  A  216  U  224 000 ----  ----   cWW  AAA   179   171     8                                                                             -     -     -
    1S72      0.4487  A 2840  C 2047  G 1728 000 ----  ----   cWW  AAA   682   989   307                                                                             -     -     -
""")

FR3D_ADJACENT = StringIO("""
Filename Discrepancy       1       2       3 Cha   1-2   1-3   2-3 Con   1-2   1-3   2-3   1-2   1-3   2-1   2-3   3-1   3-2   1-2   1-3   2-1   2-3   3-1   3-2   1-2   1-3   2-3
    1S72      0.4292  A   49  G  149  C   42 000 ----  ----   cWW  AAA    98     7   105                                                                             -     -     -
    1S72      0.4358  A   99  C   82  G   92 000 ----  ----   cWW  AAA    17     7    10                                                                             -     -     -    
""")

FR3D_NO_STEM = StringIO("""
Filename Discrepancy       1       2       3 Cha   1-2   1-3   2-3 Con   1-2   1-3   2-3   1-2   1-3   2-1   2-3   3-1   3-2   1-2   1-3   2-1   2-3   3-1   3-2   1-2   1-3   2-3
    1S72      0.4590  A 2019  C 1830  G  820 000 ----  ----   ncWW AAS   177  1158   981                                      2BR                                    -     -     -
""")
unittest.TestCase.longMessage = True
class TestAMinor(unittest.TestCase):
    def setUp(self):
        A_MULTICHAIN_FR3D_OUTPUT.seek(0)
        A_FR3D_OUTPUT.seek(0)
        FR3D_ADJACENT.seek(0)
        FR3D_NO_STEM.seek(0)
        ANOTHER_FR3D_OUTPUT.seek(0)

    def test_fr3d_orientation(self):
        cg = ftmc.CoarseGrainRNA("test/fess/data/1S72_0.cg")
        a = fba.parse_fred(30, {"1S72":cg}, ANOTHER_FR3D_OUTPUT)
        for geo1, geo2 in it.combinations(a, 2):
            self.assertLess(abs(geo1.dist - geo2.dist), 20, msg = "Different FR3D hits should have similar coarse-grained geometry: dist. {}, {}".format(geo1, geo2))
            self.assertLess(abs(geo1.angle1 - geo2.angle1), 2., msg = "Different FR3D hits should have similar coarse-grained geometry: angle1. {}, {}".format(geo1, geo2))
            self.assertLess(abs(geo1.angle2 - geo2.angle2), 2., msg = "Different FR3D hits should have similar coarse-grained geometry: angle2. {}, {}".format(geo1, geo2))
    def test_get_relative_orientation(self):
        cg = ftmc.CoarseGrainRNA("test/fess/data/1S72_0.cg")
        for loop in cg.edges["s0"]:
            d, a1, a2 = fba.get_relative_orientation(cg, kloop, "s0")
            self.assertEqual(d, 0)
            
    def test_get_relative_orientation2(self):
        cg = ftmc.CoarseGrainRNA(dotbracket_str = "...(((...)))...", seq="AAAGGGAAACCCAAA")
        cg.coords["s0"] = [0,0,0.], [0,0,10.]
        cg.twists["s0"] = [0,1,0.], [0, -1., 0]
        cg.coords["f0"] = [5.,0,0], [10.,0,0]
        cg.coords["t0"] = [5.,0,5], [10., 0, 5]
        cg.coords["h0"] = [5.,0,10],[10.,0,10]
        d, a1, a2 = fba.get_relative_orientation(cg, "f0", "s0")
        self.assertAlmostEqual(d, 5)
        self.assertAlmostEqual(a1, math.pi/2)
        self.assertAlmostEqual(a2, math.pi/2)
        d, a1, a2 = fba.get_relative_orientation(cg, "h0", "s0")
        self.assertAlmostEqual(d, 5)
        self.assertAlmostEqual(a1, math.pi/2)
        self.assertAlmostEqual(a2, math.pi/2)
        d, a1, a2 = fba.get_relative_orientation(cg, "t0", "s0")
        self.assertAlmostEqual(d, 5)
        self.assertAlmostEqual(a1, math.pi/2)
        self.assertAlmostEqual(a2, 0)

        
    def test_parse_fred_missing_chain(self):
        cg = ftmc.CoarseGrainRNA("test/fess/data/1S72_0.cg")
        a = fba.parse_fred(30, {"1S72":cg}, A_MULTICHAIN_FR3D_OUTPUT)
        #Chain 9 is not connected to chain 0, thus not present in the cg-file.
        self.assertEqual(len(a), 0)
        
    def test_parse_fred_adjacent(self):
        cg = ftmc.CoarseGrainRNA("test/fess/data/1S72_0.cg")
        a = fba.parse_fred(30, {"1S72":cg}, FR3D_ADJACENT)
        #The loop is adjacent to the stem.
        self.assertEqual(len(a), 0)
    def test_parse_fred_no_stem(self):
        cg = ftmc.CoarseGrainRNA("test/fess/data/1S72_0.cg")
        a = fba.parse_fred(30, {"1S72":cg}, FR3D_NO_STEM)
        #The receptor is no canonical stem.
        self.assertEqual(len(a), 0)
    def test_parse_fred1(self):
        cg = ftmc.CoarseGrainRNA("test/fess/data/1S72_0.cg")
        a = fba.parse_fred(30, {"1S72":cg}, A_FR3D_OUTPUT)
        self.assertEqual(len(a), 1)
        geometry, = a
        self.assertEqual(geometry.pdb_id, "1S72")
    def test_parse_fred2(self):
        cg = ftmc.CoarseGrainRNA("test/fess/data/1S72_0.cg")
        a = fba.parse_fred(30, {"1S72":cg}, ANOTHER_FR3D_OUTPUT)
        self.assertEqual(len(a), 4)
