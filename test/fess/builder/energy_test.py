from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)

import unittest, sys
import numpy as np
import numpy.testing as nptest
import forgi.threedee.model.projection2d as ftmp
import fess.builder.energy as fbe
import fess.builder.models as fbm
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.vector as ftuv
class TestProjectionMatchEnergySetup(unittest.TestCase):
    def test_ProjectionMatchEnergy_init(self):
        try:
            energy=fbe.ProjectionMatchEnergy({("h1","h2"):15})
        except Exception as e: 
            assert False, "Error during init of projectionMatchEnergy, {}".format(e)

class TestProjectionMatchEnergy(unittest.TestCase):
    def setUp(self):
        cg1 = ftmc.CoarseGrainRNA('test/fess/data/1GID_A-structure1.coord')
        cg2 = ftmc.CoarseGrainRNA('test/fess/data/1GID_A-structure2.coord')  
        self.sm1 = fbm.SpatialModel(cg1)
        self.sm2 = fbm.SpatialModel(cg2)            
        self.sm1.load_sampled_elems()
        self.sm2.load_sampled_elems()
        self.energy1a=fbe.ProjectionMatchEnergy({("h0","m0"):23.26, ("h0","h1"):43.27, ("h1","m0"):38.13})
        self.energy1b=fbe.ProjectionMatchEnergy({("h0","m2"):37.2, ("h0","i4"):40.14, ("h0","h1"):40.14, ("m2","i4"):3.0,("m2","h1"):6.75, ("i4","h1"):6.91})

        self.energy2a=fbe.ProjectionMatchEnergy({("h1","h2"):42.74, ("h1","s9"):56.40, ("h1","m0"):51.96, ("h2","m0"):19.98,("h2","s9"):18.62, ("m0","s9"):7.85})
        self.energy2b=fbe.ProjectionMatchEnergy({("h2","h1"):47.95, ("h2","m3"):34.14, ("h2","i1"):21.67, ("h1","m3"):13.82,("i1","h1"):29.33, ("i1","m3"):16.52})

        self.energyA1=fbe.ProjectionMatchEnergy({("h0","h1"):40.49, ("h0", "m0"): 39.78, ("h0","t1"):43.70, ("h0", "i6"):32.64, ("h1","m0"):51.15, ("h1","t1"):50.02,("h1", "i6"):54.61, ("m0","t1"):6.08, ("m0","i6"):12.75, ("t1","i6"):18.82})
        self.energyA2=fbe.ProjectionMatchEnergy({("h0","h1"):40.49, ("h0", "m0"): 39.78, ("h0","t1"):43.70, ("h1","m0"):51.15, ("h1","t1"):50.02, ("m0","t1"):6.08})
        self.energyA3=fbe.ProjectionMatchEnergy({("h0","h1"):40.49, ("h0", "m0"): 39.78, ("h1","m0"):51.15})
        self.energyA4=fbe.ProjectionMatchEnergy({("h0", "m0"): 39.78})
        return
    def test_ProjectionMatchEnergy_eval_energy_correct_projection(self):
        ENERGY_TOLERANCE=0.2
        VECTOR_A_TOLERANCE=0.05
        e=self.energy1a.eval_energy(self.sm1)
        self.assertLessEqual(e, ENERGY_TOLERANCE)
        targetdir=np.array([0.362,0.023, -0.826])
        targetdir=targetdir/ftuv.magnitude(targetdir)
        if self.energy1a.projDir[2]>0:
            targetdir=-1*targetdir
        nptest.assert_allclose(self.energy1a.projDir, targetdir,atol=VECTOR_A_TOLERANCE)
        e=self.energy1b.eval_energy(self.sm1)
        self.assertLessEqual(e, ENERGY_TOLERANCE)
        targetdir= np.array([-0.193,-0.319,0.074])
        targetdir=targetdir/ftuv.magnitude(targetdir)
        if self.energy1b.projDir[1]>0:
            targetdir=-1*targetdir
        nptest.assert_allclose(self.energy1b.projDir, targetdir, atol=VECTOR_A_TOLERANCE)
        e=self.energy2a.eval_energy(self.sm2)
        self.assertLessEqual(e, ENERGY_TOLERANCE)
        targetdir=np.array([-0.223,0.048,-0.579])
        targetdir=targetdir/ftuv.magnitude(targetdir)
        if self.energy2a.projDir[2]>0:
            targetdir=-1*targetdir
        nptest.assert_allclose(self.energy2a.projDir, targetdir,atol=VECTOR_A_TOLERANCE)
        e=self.energy2b.eval_energy(self.sm2)
        self.assertLessEqual(e, ENERGY_TOLERANCE)
        targetdir=np.array([-0.464,-0.345,-0.192])
        targetdir=targetdir/ftuv.magnitude(targetdir)
        if self.energy2b.projDir[2]>0:
            targetdir=-1*targetdir
        nptest.assert_allclose(self.energy2b.projDir, targetdir, atol=VECTOR_A_TOLERANCE)

    def test_ProjectionMatchEnergy_eval_energy_wrong_projection(self):  
        WRONG_ENERGY=1.8
        e=self.energy1a.eval_energy(self.sm2)
        self.assertGreater(e,WRONG_ENERGY)
        e=self.energy1b.eval_energy(self.sm2)
        self.assertGreater(e,WRONG_ENERGY)
        e=self.energy2a.eval_energy(self.sm1)
        self.assertGreater(e,WRONG_ENERGY)
        e=self.energy2b.eval_energy(self.sm1)
        self.assertGreater(e,WRONG_ENERGY)

    def test_ProjectionMatchEnergy_eval_energy_effectOf_num_constraints(self):
        """
        Adding more constraints should not increase the energy too much

        More constraints mean that some constraints are likely to be not fulfilled, which can lead to an increase of energy.
        However, if they are equally fulfilled, the energy should stay approximately the same.
        """
        ENERGY_CHANGE=1.2
        #Correct 3D structure
        e1=self.energyA1.eval_energy(self.sm2)
        e2=self.energyA2.eval_energy(self.sm2)
        e3=self.energyA3.eval_energy(self.sm2)
        e4=self.energyA4.eval_energy(self.sm2)
        #print("ENERGIES", e1,e2,e3,e4)
        self.assertLess(abs(e2-e1)/e1,ENERGY_CHANGE)
        self.assertLess(abs(e3-e1)/e1,ENERGY_CHANGE)
        self.assertLess(abs(e3-e2)/e2,ENERGY_CHANGE)
        self.assertLess(abs(e4-e1)/e1,ENERGY_CHANGE)
        self.assertLess(abs(e4-e2)/e2,ENERGY_CHANGE)
        self.assertLess(abs(e4-e3)/e3,ENERGY_CHANGE)
        #"Wrong" 3D structure
        e1=self.energyA1.eval_energy(self.sm1)
        e2=self.energyA2.eval_energy(self.sm1)
        e3=self.energyA3.eval_energy(self.sm1)
        e4=self.energyA4.eval_energy(self.sm1)
        #print("ENERGIES", e1,e2,e3,e4)
        self.assertLess(abs(e3-e2)/e2,ENERGY_CHANGE)
        self.assertLess(abs(e2-e1)/e1,ENERGY_CHANGE)
        self.assertLess(abs(e3-e1)/e1,ENERGY_CHANGE)
        self.assertLess(abs(e4-e1)/e1,4)
        self.assertLess(abs(e4-e2)/e2,3)
        self.assertLess(abs(e4-e3)/e3,ENERGY_CHANGE)


