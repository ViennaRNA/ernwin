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

try:
    from unittest.mock import Mock #python3
except:
    from mock import Mock
    
def add_stem_coordinates(cg, stem, start, direction=[0.,0.,10.]):
    """
    :param cg: The CoarseGrainedRNA
    :param stem: STRING, e.g. "s0"
    :param start: A sequence of 3 floats (coordinates)
    :param direction: A vector pointing from the start to the end
    """
    cg.coords[stem] = start, np.array(start)+np.array(direction)
    cg.twists[stem] = (ftuv.get_orthogonal_unit_vector(cg.coords.get_direction(stem)), 
                      -ftuv.get_orthogonal_unit_vector(cg.coords.get_direction(stem)))

class TestClashEnergy(unittest.TestCase):
    def setUp(self):
        self.cg=ftmc.CoarseGrainRNA('test/fess/data/1GID_A-structure1.coord')
        self.cg2=ftmc.CoarseGrainRNA('test/fess/data/1GID_A-structure2.coord')
        self.cg_clash=ftmc.CoarseGrainRNA('test/fess/data/1GID_A-clash.coord')
        self.cg.add_all_virtual_residues()
        self.cg2.add_all_virtual_residues()
        self.cg_clash.add_all_virtual_residues()
        #self.sm = fbm.SpatialModel(cg)
        #self.sm.load_sampled_elems()
        #self.sm2 = fbm.SpatialModel(cg2)
        #self.sm2.load_sampled_elems()
        #self.sm_clash= fbm.SpatialModel(cg_clash)
        #self.sm_clash.load_sampled_elems()
        self.energy=fbe.StemVirtualResClashEnergy()
    def test_stem_virtual_res_clash_energy_with_nodes(self):
        self.assertEqual(self.energy.eval_energy(self.cg), 0.)
        nodes=['h2', 'h0', 'h1', 's9', 's8', 's3', 's2', 's1', 's0', 's7', 's6', 's5', 's4', 'm1', 
               'm0', 'm3', 'm2', 's11', 's10', 'i1', 'i0', 'i3', 'i2', 'i5', 'i4', 'i7', 'i6', 't1']
        self.assertEqual(self.energy.eval_energy(self.cg, nodes=nodes), 0.)
        nodes=['s9', 's8', 's3', 's2', 's1', 's0', 's7', 's6', 's5', 's4', 's11', 's10']
        self.assertEqual(self.energy.eval_energy(self.cg, nodes=nodes), 0.)
        #Raises, if not at least one stem in Nodes.
        with self.assertRaises(ValueError):
            self.assertEqual(self.energy.eval_energy(self.cg, nodes=["i5"]), 0.)
        with self.assertRaises(ValueError):
            self.assertEqual(self.energy.eval_energy(self.cg, nodes=[]), 0.)
        #Raises, if stem not in Graph
        with self.assertRaises(KeyError):
            self.assertEqual(self.energy.eval_energy(self.cg, nodes=["s220"]), 0.)
        #Structure with a clash
        self.assertGreater(self.energy.eval_energy(self.cg_clash), 100.)
        self.assertGreater(self.energy.eval_energy(self.cg_clash, nodes=["s7", "s11"]), 100.)

        
    def test_energy_independent_of_nodes(self):
        for i,cg in enumerate([self.cg, self.cg_clash, self.cg2]):                
            e=self.energy.eval_energy(cg)
            for l in range(10,len(cg.defines),2):
                nodes=random.sample(cg.defines.keys(),l)
                try:
                    e_nodes=self.energy.eval_energy(cg, nodes=nodes)
                except ValueError: #No stem in nodes
                    e_nodes=0
                np.set_printoptions(threshold=np.nan)
                self.assertLessEqual(e_nodes, e, "{} is not <= {}. The clash energy should be "
                                     "smaller or the same, if nodes are used. Nodes used were {} "
                                     "for spatial model {}.".format(e_nodes, e, nodes, i))

    def test_bad_bulges(self):
        self.energy.eval_energy(self.cg_clash)
        self.assertEqual(set(self.energy.bad_bulges), {"s7", "s11"})

class TestJunctionConstraintEnergy(unittest.TestCase):
    def setUp(self):
        self.cg = ftmc.CoarseGrainRNA('test/fess/data/1GID_A.cg')
        self.cg.add_all_virtual_residues()
        self.junction_energy = fbe.RoughJunctionClosureEnergy()
        self.cg_bad = ftmc.CoarseGrainRNA(dotbracket_str = '(((...)))...(((...)))', seq = "AAAGGGUUUGGGUUUGGGAAA")
        add_stem_coordinates(self.cg_bad, "s0", [0.,0.,0.], [0.,0.,10.])
        self.cg_bad.coords["h0"] = self.cg_bad.coords["s0"][1], self.cg_bad.coords["s0"][1]+[3.,6.,0.]
        self.cg_bad.coords["m0"] = self.cg_bad.coords["s0"][0], self.cg_bad.coords["s0"][0]+[2.,1.,300.]
        add_stem_coordinates(self.cg_bad, "s1", self.cg_bad.coords["m0"][1], [12.,2.,-2.])
        self.cg_bad.coords["h1"] = self.cg_bad.coords["s1"][1], self.cg_bad.coords["s1"][1]+[1.,6.,0.]
        self.cg_bad.add_all_virtual_residues()
    def test_junction_energy_ok(self):
        self.assertEqual(self.junction_energy.eval_energy(self.cg), 0)
    def test_junction_energy_bad(self):
        self.assertGreater(self.junction_energy.eval_energy(self.cg_bad), 1000)
        self.assertIn("m0", self.junction_energy.bad_bulges)    
    def test_junction_energy_nodes(self):
        self.assertEqual(self.junction_energy.eval_energy(self.cg, nodes=["s0", "h0", "s1"]), 0)
        self.assertGreater(self.junction_energy.eval_energy(self.cg_bad, nodes=["m0"]), 1000)
    def test_energy_independent_of_nodes(self):
        cg = ftmc.CoarseGrainRNA('test/fess/data/4way.cg')
        cg.add_all_virtual_residues()
        self.assertEqual(self.junction_energy.eval_energy(cg, nodes=["m0"]), self.junction_energy.eval_energy(cg))
class TestSLDEnergies(unittest.TestCase):
    def setUp(self):
        self.cg1 = ftmc.CoarseGrainRNA('test/fess/data/1GID_A.cg')
        self.cg1.add_all_virtual_residues()
        self.cg_far = ftmc.CoarseGrainRNA(dotbracket_str = '(((...)))...(((...)))', seq = "AAAGGGUUUGGGUUUGGGAAA")
        add_stem_coordinates(self.cg_far, "s0", [0.,0.,0.], [0.,0.,10.])
        self.cg_far.coords["h0"] = self.cg_far.coords["s0"][1], self.cg_far.coords["s0"][1]+[3.,6.,0.]
        self.cg_far.coords["m0"] = self.cg_far.coords["s0"][0], self.cg_far.coords["s0"][0]+[2.,1.,200.]
        add_stem_coordinates(self.cg_far, "s1", self.cg_far.coords["m0"][1], [12.,2.,-2.])
        self.cg_far.coords["h1"] = self.cg_far.coords["s1"][1], self.cg_far.coords["s1"][1]+[1.,6.,0.]
        self.cg_far.add_all_virtual_residues()
    
        self.cg_five = ftmc.CoarseGrainRNA(dotbracket_str = '(((...)))...(((...)))...(((...)))...(((...)))...(((...)))', seq = "AAAGGGUUUGGGUUUGGGAAAGGGUUUGGGAAAGGGUUUGGGAAAGGGUUUGGGAAA")
        add_stem_coordinates(self.cg_five, "s0", [0.,0.,0.], [0.,0.,10.])
        for i in range(1, 5):
            m="m{}".format(i-1)
            prev_s = "s{}".format(i-1)
            s = "s{}".format(i)
            self.cg_five.coords[m] = self.cg_five.coords[prev_s][0], self.cg_five.coords[prev_s][0]+[0.,3.+i,0.]
            add_stem_coordinates(self.cg_five, s, self.cg_five.coords[m][1], [0.,0.,10.])        
        for i in range(5):
            s="s{}".format(i)
            h="h{}".format(i)
            self.cg_five.coords[h] = self.cg_five.coords[s][1], self.cg_five.coords[s][1]+[3.,0.,0.]
        self.cg_five.add_all_virtual_residues()

    
    
    
    def test_SDL_1GID(self):
        energy = fbe.ShortestLoopDistancePerLoop.from_cg(self.cg1, None, None)
        self.assertGreater(energy.eval_energy(self.cg1, background=True), -30)
        self.assertLess(energy.eval_energy(self.cg1, background=True), 30)
        self.assertGreater(energy.eval_energy(self.cg1, background=False), 0)
        self.assertLess(energy.eval_energy(self.cg1, background=False), 50)
        
    def test_SLD_far(self):
        energy = fbe.ShortestLoopDistancePerLoop.from_cg(self.cg_far, None, None)
        # At first, we do not use any artificial data
        energy._lsp_data_weight = 1
        energy._lsp_artificial_weight = 0
        energy.kde_resampling_frequency = 1
        energy.reset_kdes(self.cg_far.seq_length)
        
        e1 = energy.eval_energy(self.cg_far)
        self.assertLess(e1, -500)
        energy.accept_last_measure()
        
        e2 = energy.eval_energy(self.cg_far)
        self.assertGreater(e2, 1000)
        energy.accept_last_measure()
        
        e3 = energy.eval_energy(self.cg_far)
        self.assertGreater(e3, 1000)
        
        
        # Now we use artificial data. Energy values should be less extreme
        energy._lsp_data_weight = 3
        energy._lsp_artificial_weight = 1
        energy.reset_kdes(self.cg_far.seq_length)
        
        e1 = energy.eval_energy(self.cg_far)#, plot_debug = True)
        self.assertLess(e1, 100)
        self.assertGreater(e1, -100)
        energy.accept_last_measure()
        
        e2 = energy.eval_energy(self.cg_far)
        self.assertLess(e2, 100)
        self.assertGreater(e2, -100)
        self.assertGreater(e2, e1)
        energy.accept_last_measure()
        
        e3 = energy.eval_energy(self.cg_far)
        self.assertLess(e3, 100)
        self.assertGreater(e3, -100)

    def test_minimal_h_h_distance(self):
        self.assertEqual(fbe._minimal_h_h_distance(self.cg_five, "h0", self.cg_five.hloop_iterator()), 4.)
        self.assertEqual(fbe._minimal_h_h_distance(self.cg_five, "h1", self.cg_five.hloop_iterator()), 4.)
        self.assertEqual(fbe._minimal_h_h_distance(self.cg_five, "h2", self.cg_five.hloop_iterator()), 5.)
        self.assertEqual(fbe._minimal_h_h_distance(self.cg_five, "h3", self.cg_five.hloop_iterator()), 6.)
        self.assertEqual(fbe._minimal_h_h_distance(self.cg_five, "h4", self.cg_five.hloop_iterator()), 7.)
        self.assertEqual(fbe._minimal_h_h_distance(self.cg_five, "h0", ["h3", "h4"]), 15.)

class TestAMinorEnergy(unittest.TestCase):
    def setUp(self):
        self.cg = ftmc.CoarseGrainRNA('test/fess/data/1GID_A.cg')
        self.cg.add_all_virtual_residues()
        self.energy = fbe.AMinorEnergy.from_cg(self.cg, None, None)
    def test_AME_energy(self):
        e = self.energy.eval_energy(self.cg)#, plot_debug = True)
        self.assertLess(e, 0)
        self.assertGreater(e, -200)

class TestABCs(unittest.TestCase):
    def setUp(self):
        self.cg = ftmc.CoarseGrainRNA('test/fess/data/1GID_A.cg')
        self.cg.add_all_virtual_residues()
    def test_background_kde_is_resampled(self):
        e = fbe.RadiusOfGyrationEnergy(self.cg.seq_length, prefactor=30)
        e._resample_background_kde = Mock()
        e.accept_last_measure()
        e._resample_background_kde.assert_not_called()
        e.accept_last_measure()
        e._resample_background_kde.assert_not_called()
        e.accept_last_measure()
        e._resample_background_kde.assert_called_once_with()
        e._resample_background_kde.reset_mock()
        e.accept_last_measure()
        e._resample_background_kde.assert_not_called()
        e._resample_background_kde.reset_mock()
        e.kde_resampling_frequency = 1
        e.accept_last_measure()
        e._resample_background_kde.assert_called_once_with()


class TestCombinedEnergy(unittest.TestCase):
    def test_getattr(self):
        e = fbe.CombinedEnergy()
        e.accept_last_measure() #Does nothing
        e.energies.append(fbe.NormalDistributedRogEnergy(12, 35))
        e.accept_last_measure() #Calls accept_last_measure fof NDR energy
        mock_energy = Mock()
        mock_energy.accept_last_measure = Mock()
        e.energies.append(mock_energy)
        e.accept_last_measure() #Calls accept_last_measure of NDR and mock_energy
        mock_energy.accept_last_measure.assert_called_once_with()
        
        with self.assertRaises(AttributeError):
            e.do_something_else()
    def test_hasinstance(self):
        e = fbe.CombinedEnergy([1, 1.2]) #Thanks to Ducktyping, energies can be any object if we don't use them
        self.assertTrue(e.hasinstance(int))
        self.assertTrue(e.hasinstance(float))
        self.assertFalse(e.hasinstance(str))
    def test_hasinstance_nested(self):
        e = fbe.CombinedEnergy([fbe.CombinedEnergy([1, 1.2])]) #Thanks to Ducktyping, energies can be any object if we don't use them
        self.assertTrue(e.hasinstance(int))
        self.assertTrue(e.hasinstance(float))
        self.assertFalse(e.hasinstance(str))


            
class TestGyrationRadiusEnergies(unittest.TestCase):
    def setUp(self):
        self.cg = ftmc.CoarseGrainRNA('test/fess/data/1GID_A.cg')
        self.cg.add_all_virtual_residues()

    def test_ROG_energy(self):
        energyfunction = fbe.RadiusOfGyrationEnergy(self.cg.seq_length, prefactor=30)
        energy = energyfunction.eval_energy(self.cg, background = False)
        self.assertLess(energy, 100)
        self.assertGreater(energy, 0)
        energy = energyfunction.eval_energy(self.cg, background = True)
        self.assertLess(energy, 50)
        self.assertGreater(energy, -50)
        
    def test_NDR_energy(self):
        energyfunction = fbe.NormalDistributedRogEnergy(self.cg.seq_length, 35)
        energyBG = energyfunction.eval_energy(self.cg, background = True)
        self.assertLess(energyBG, 1000)
        self.assertGreater(energyBG, -1000)

    def test_ROG_energy_last_measure(self):
        energyfunction = fbe.RadiusOfGyrationEnergy(self.cg.seq_length)
        energy = energyfunction.eval_energy(self.cg)
        energyfunction.accept_last_measure()
        self.assertEqual(energyfunction.accepted_measures[-1], self.cg.radius_of_gyration("vres"))
        
class TestProjectionMatchEnergySetup(unittest.TestCase):
    def test_ProjectionMatchEnergy_init(self):
        try:
            energy=fbe.ProjectionMatchEnergy({("h1","h2"):15})
        except Exception as e: 
            assert False, "Error during init of projectionMatchEnergy, {}".format(e)
            
@unittest.skip("Projection match energy: The 3D structures changed, so we need to update the tests.")
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


class TestConvenienceFunctions(unittest.TestCase):
    def test_energies_from_string_single(self):
        cg = ftmc.CoarseGrainRNA('test/fess/data/1GID_A-structure1.coord')
        
        energy = fbe.energies_from_string("SLD", cg)
        self.assertTrue(energy.hasinstance(fbe.ShortestLoopDistancePerLoop))
    
    def test_from_cg_AME(self):
        # A in hairpin and interiorloop
        cg1 = ftmc.CoarseGrainRNA(dotbracket_str = '(((..(((...)))..)))', seq = "GGGAAGGGAAACCCAACCC")
        e1 = fbe.AMinorEnergy.from_cg(cg1, None, None)
        self.assertEqual(len(e1.energies), 2)
        # A only in hairpin
        cg2 =  ftmc.CoarseGrainRNA(dotbracket_str = '(((..(((...)))..)))', seq = "GGGUUGGGAAACCCUUCCC")
        e2 = fbe.AMinorEnergy.from_cg(cg2, None, None)
        self.assertEqual(len(e2.energies), 1)
        self.assertEqual(e2.energies[0].loop_type, 0)
        # A only in IL
        cg2 =  ftmc.CoarseGrainRNA(dotbracket_str = '(((..(((...)))..)))', seq = "GGGUAGGGUUUCCCUUCCC")
        e2 = fbe.AMinorEnergy.from_cg(cg2, None, None)
        self.assertEqual(len(e2.energies), 1)
        self.assertEqual(e2.energies[0].loop_type, 1)
        
    def test_energies_from_string_multiple(self):
        cg = ftmc.CoarseGrainRNA('test/fess/data/1GID_A-structure1.coord')

        energy = fbe.energies_from_string("AME,ROG", cg)
        self.assertTrue(energy.hasinstance(fbe.AMinorEnergy))
        self.assertTrue(energy.hasinstance(fbe.RadiusOfGyrationEnergy))
        
    def test_energies_from_string_prefactor(self):
        cg = ftmc.CoarseGrainRNA('test/fess/data/1GID_A-structure1.coord')

        energy = fbe.energies_from_string("20AME,15NDR,ROG", cg)
        for e in energy.iterate_energies():
            if isinstance(e, fbe.AMinorEnergy):
                self.assertEqual(e.prefactor, 20)
            elif isinstance(e, fbe.NormalDistributedRogEnergy):
                self.assertEqual(e.prefactor, 15)
            elif isinstance(e, fbe.RadiusOfGyrationEnergy):
                self.assertEqual(e.prefactor, fbe.DEFAULT_ENERGY_PREFACTOR)
    def test_energies_from_string_adjustment(self):
        cg = ftmc.CoarseGrainRNA('test/fess/data/1GID_A-structure1.coord')

        energy = fbe.energies_from_string("AME1.5,NDR,ROG3.2", cg)
        for e in energy.iterate_energies():
            if isinstance(e, fbe.AMinorEnergy):
                self.assertEqual(e.adjustment, 1.5)
            elif isinstance(e, fbe.NormalDistributedRogEnergy):
                self.assertEqual(e.adjustment, 1.)
            elif isinstance(e, fbe.RadiusOfGyrationEnergy):
                self.assertEqual(e.adjustment, 3.2)
        
    def test_energies_from_string_simulated_annealing_2(self):
        cg = ftmc.CoarseGrainRNA('test/fess/data/1GID_A-structure1.coord')
        energy = fbe.energies_from_string("AME1.5_2.4,NDR,20_39ROG", cg, 1000)
        for e in energy.iterate_energies():
            if isinstance(e, fbe.AMinorEnergy):
                self.assertEqual(e.adjustment, 1.5)
                self.assertEqual(e._adj_stepwidth, 0.1)
                self.assertEqual(e._adj_update_freq, 100)
            elif isinstance(e, fbe.NormalDistributedRogEnergy):
                self.assertEqual(e._adj_update_freq, float("inf"))
            elif isinstance(e, fbe.RadiusOfGyrationEnergy):
                self.assertEqual(e.prefactor, 20)
                self.assertEqual(e._pf_stepwidth, 1)
                self.assertEqual(e._pf_update_freq, 50)
                
    def test_energies_from_string_simulated_annealing_3(self):
        cg = ftmc.CoarseGrainRNA('test/fess/data/1GID_A-structure1.coord')
        energy = fbe.energies_from_string("AME1.5_0.2_2.3,NDR,20_5_35ROG", cg, 1000)
        for e in energy.iterate_energies():
            if isinstance(e, fbe.AMinorEnergy):
                self.assertEqual(e.adjustment, 1.5)
                self.assertEqual(e._adj_stepwidth, 0.2)
                self.assertEqual(e._adj_update_freq, 200)
            elif isinstance(e, fbe.NormalDistributedRogEnergy):
                self.assertEqual(e._adj_update_freq, float("inf"))
            elif isinstance(e, fbe.RadiusOfGyrationEnergy):
                self.assertEqual(e.prefactor, 20)
                self.assertEqual(e._pf_stepwidth, 5)
                self.assertEqual(e._pf_update_freq, 250)
                
    def test_energies_from_string_kwargs_clamp(self):
        cg = ftmc.CoarseGrainRNA('test/fess/data/1GID_A-structure1.coord')
        energy = fbe.energies_from_string("CLA", cg, 1000, cla_pairs = [("s0", "s1"),("s0", "h2")])
        for e in energy.iterate_energies():
            self.assertEqual(e.adjustment, e._DEFAULT_CLAMP_DISTANCE)
            self.assertEqual(e.from_elem, "s0")
            self.assertIn(e.to_elem, ["s1", "h2"])
        with self.assertRaises(TypeError) as ex: #Missing kwarg for CLA energy
            energy = fbe.energies_from_string("CLA", cg, 1000)
        self.assertIn("CLA", str(ex.exception))
        self.assertIn("cla_pairs", str(ex.exception))

class TestHelperFunctions(unittest.TestCase):
    def test__iter_subgraphs(self):
        cg = ftmc.CoarseGrainRNA(dotbracket_str = "...(((...)))...(((...)))...(((...(((...)))...)))", seq="AAAGGGAAACCCAAAGGGAAACCCAAAGGGUUUGGGAAACCCUUUCCC")
        sgs = fbe._iter_subgraphs(cg, 1)
        self.assertEqual(len(sgs), 1)
        sgs = fbe._iter_subgraphs(cg, True)
        self.assertGreater(len(sgs), 4)
