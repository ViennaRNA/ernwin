import os.path as op

import unittest

import fess.builder.energy as fbe

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.stats as ftms
import forgi.utilities.debug as fud

import fess.builder.models as fbm

class TestModel(unittest.TestCase):
    def setUp(self):
        self.cg = ftmc.CoarseGrainRNA('test/fess/data/1GID_A.cg')
        self.real_stats_fn = 'test/fess/data/real.stats'
        self.filtered_stats_fn = 'test/fess/data/filtered_stats_1gid.csv'

        self.conf_stats = ftms.ConformationStats(self.real_stats_fn)
        self.filtered_stats = ftms.FilteredConformationStats(self.real_stats_fn, 
                                                             self.filtered_stats_fn)
        
        self.sm = fbm.SpatialModel(self.cg, conf_stats = self.conf_stats)
        self.sm.sample_stats()

        return

    def test_filtered_traverse_and_build(self):
        fcs = self.filtered_stats
        ftms.set_conformation_stats(fcs)
        sm = fbm.SpatialModel(self.cg, conf_stats=fcs)
        sm.sample_stats()

        sm.traverse_and_build()
        sm.bg.to_file('temp.cg')

    def test_sample_elems(self):
        sm = fbm.SpatialModel(self.cg, conf_stats = self.conf_stats)

        sm.sample_stats()
    
    def test_get_random_stem_stats(self):
        self.sm.get_random_stem_stats('s0')

    def test_traverse_and_build(self):
        return
        sm = fbm.SpatialModel(self.cg, conf_stats=self.conf_stats)
        sm.sample_stats()

        sm.traverse_and_build()

        cg = ftmc.CoarseGrainRNA('test/fess/data/1ymo_pk.cg')
        sm = fbm.SpatialModel(cg, conf_stats=self.conf_stats)
        sm.sample_stats()
        sm.traverse_and_build()
        sm.bg.to_file('temp1.cg')

        #pseudoknot
        cg = ftmc.CoarseGrainRNA(op.expanduser('~/doarse/4LVV_A/temp.cg'))
        cg = ftmc.from_pdb(op.expanduser('~/doarse/4LVV_A/temp.pdb'),
                           remove_pseudoknots=False)
        sm = fbm.SpatialModel(cg, conf_stats=self.conf_stats)
        sm.sample_stats()
        sm.traverse_and_build()
        #sm.bg.to_file('temp1.cg')

    def test_new_traverse_and_build(self):
        import time, sys

        time1 = time.time()
        for i in range(10):
            cg = ftmc.CoarseGrainRNA('test/fess/data/4P8Z.cg')
            sm = fbm.SpatialModel(cg)
            sm.sample_stats()

            sm.constraint_energy = fbe.CombinedEnergy([fbe.RoughJunctionClosureEnergy(), 
                                                    fbe.CoarseStemClashEnergy(), 
                                                    fbe.StemVirtualResClashEnergy()])
            sm.traverse_and_build()
        time2 = time.time()
        print >>sys.stderr, "traverse_and_build, elapsed_time:", time2 - time1

        time1 = time.time()
        for i in range(10):
            cg = ftmc.CoarseGrainRNA('test/fess/data/4P8Z.cg')
            sm = fbm.SpatialModel(cg)
            sm.sample_stats()

            sm.constraint_energy = fbe.CombinedEnergy([fbe.RoughJunctionClosureEnergy(), 
                                                    fbe.CoarseStemClashEnergy(), 
                                                    fbe.StemVirtualResClashEnergy()])
            sm.new_traverse_and_build()
        time2 = time.time()

        print >>sys.stderr, "new_traverse_and_build, elapsed time:", time2 - time1
