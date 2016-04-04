from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)

import unittest, sys, random
import numpy as np
import numpy.testing as nptest
import forgi.projection.projection2d as ftmp
import fess.builder.energy as fbe
import fess.builder.models as fbm
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.vector as ftuv

class TestClashEnergy(unittest.TestCase):
    def setUp(self):
        cg=ftmc.CoarseGrainRNA('test/fess/data/someSTMV.coord')
        cg2=ftmc.CoarseGrainRNA('test/fess/data/anotherSTMV.coord')
        cg_clash=ftmc.CoarseGrainRNA('test/fess/data/1GID_A-clash.coord')
        self.sm = fbm.SpatialModel(cg)
        self.sm.load_sampled_elems()
        self.sm2 = fbm.SpatialModel(cg2)
        self.sm2.load_sampled_elems()
        self.sm_clash= fbm.SpatialModel(cg_clash)
        self.sm_clash.load_sampled_elems()
        self.energy=fbe.StemVirtualResClashEnergy()
    def test_stem_virtual_res_clash_energy(self):
        self.assertEqual(self.energy.eval_energy(self.sm), 0.)
        nodes=['i0', 'i1', 'i10', 'i11', 'i12', 'i13', 'i14', 'i15', 'i16', 'i17', 'i18', 'i19', 
               'i2', 'i20', 'i21', 'i22', 'i23', 'i24', 'i25', 'i26', 'i27', 'i28', 'i29', 'i3', 
               'i30', 'i31', 'i32', 'i33', 'i34', 'i35', 'i36', 'i37', 'i39', 'i4', 'i40', 'i41',
               'i42', 'i5', 'i6', 'i7', 'i8', 'i9', 'm0', 'm1', 'm10', 'm11', 'm13', 'm14', 'm15',
               'm17', 'm18', 'm19', 'm20', 'm22', 'm23', 'm24', 'm25', 'm26', 'm27', 'm28', 'm29',
               'm3', 'm32', 'm33', 'm34', 'm4', 'm5', 'm6', 'm7', 'm8', 's0', 's1', 's10', 's11',
               's12', 's13', 's14', 's15', 's16', 's17', 's18', 's19', 's2', 's20', 's21', 's22',
               's23', 's24', 's25', 's26', 's27', 's28', 's29', 's3', 's30', 's31', 's32', 's33',
               's34', 's35', 's36', 's37', 's38', 's39', 's4', 's40', 's41', 's42', 's43', 's44',
               's45', 's46', 's47', 's48', 's49', 's5', 's50', 's51', 's52', 's53', 's54', 's55',
               's56', 's57', 's58', 's59', 's6', 's60', 's61', 's62', 's64', 's67', 's68', 's69', 
               's7', 's70', 's71', 's72', 's73', 's8', 's9']
        self.assertEqual(self.energy.eval_energy(self.sm, nodes=nodes), 0.)
        nodes=['s12', 's13', 's14', 's15', 's16', 's17', 's18', 's19', 's2', 's20', 's21']
        self.assertEqual(self.energy.eval_energy(self.sm, nodes=nodes), 0.)
        #Raises, if not at least one stem in Nodes.
        with self.assertRaises(ValueError):
            self.assertEqual(self.energy.eval_energy(self.sm, nodes=["i5"]), 0.)
        with self.assertRaises(ValueError):
            self.assertEqual(self.energy.eval_energy(self.sm, nodes=[]), 0.)
        #Raises, if stem not in Graph
        with self.assertRaises(KeyError):
            self.assertEqual(self.energy.eval_energy(self.sm, nodes=["s220"]), 0.)
        #Structure with a clash
        self.assertGreater(self.energy.eval_energy(self.sm_clash), 100.)
    def test_energy_independent_of_nodes(self):
        for i,sm in enumerate([self.sm, self.sm_clash, self.sm2]):
            for l in range(10,len(sm.bg.defines),2):
                nodes=random.sample(sm.bg.defines.keys(),l)
                try:
                    e_nodes=self.energy.eval_energy(sm, nodes=nodes)
                except ValueError: #No stem in nodes
                    e_nodes=0
                c1=self.energy.c
                e=self.energy.eval_energy(sm)
                c2=self.energy.c
                np.set_printoptions(threshold=np.nan)
                self.assertLessEqual(e_nodes, e, "{} is not <= {}. The clash energy should be "
                                     "smaller or the same, if nodes are used. Nodes used were {} "
                                     "for spatial model {}. VRAS1 {}, VRAS2{}".format(e_nodes, e, nodes, i, c1, c2))

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


