from __future__ import print_function, division
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, pow, range, round,
                      str, super, zip)

import os.path as op
import unittest, copy, warnings, sys
from StringIO import StringIO

#from pprint import pprint
import fess.builder.energy as fbe

import re

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.similarity as ftmsim
import forgi.threedee.model.stats as ftms
import forgi.threedee.utilities.vector as ftuv
import forgi.utilities.debug as fud

import random

import numpy as np
import fess.scripts.ernwin_new as ernwin


from mock import mock_open, patch #In python 3 this will be from unittest.mock

class AnyStringWith(str):
    def __eq__(self, other):
        return self in other
class AnyStringStartingWith(str):
    def __eq__(self, other):
        return other.startswith(self)

class WrongAssumptionAboutTestInput(ValueError):
    pass

class ErnwinTestsMixin(object):
    def assertOneWritePerIteration(self, mocked_write, iterations):
        self.longMessage = True
        num_writes = mocked_write.call_count
        comment_count = 0
        newline_count = 0
        for call in mocked_write.mock_calls:
            name, args, kwargs = call
            if args[0].startswith("#"): comment_count+=1
            if args[0] == "\n": newline_count+=1
        self.assertEqual(num_writes, iterations + 1 + comment_count + newline_count, 
                         "{} writes were recorded (ignoring {} comments, {} newlines and 1 header line). "
                         "Expected {}".format( num_writes - comment_count - 1 - newline_count,
                                               comment_count, newline_count, iterations))
    def runErnwin(self, command):

        self.exitCode = None
        open_stats = mock_open()
        open_main = mock_open()
        try:     
            with patch('sys.stdout', new_callable=StringIO) as mock_stdout:
                with patch('sys.stderr', new_callable=StringIO) as mock_stderr:
                    with patch('fess.builder.samplingStatisticsNew2.open', open_stats, create=True):
                        with patch('fess.scripts.ernwin_new.open', open_main, create=True):
                            try:
                                ernwin.main(self.parser.parse_args(command.split()[2:]))
                            except SystemExit as e:
                                self.exitCode = e.code
            if self.exitCode:
                print("ERNWIN exited with non-zero exit code {}".format(self.exitCode), file=sys.stderr)
                print("STDERR WAS:", file=sys.stderr)
                print(mock_stderr.getvalue(), file=sys.stderr)
            return open_main, open_stats, mock_stdout, mock_stderr
        except BaseException as e:
            print('Exception {}:"{}" caught'.format(type(e), e), file=sys.stderr)
            print("STDERR WAS:", file=sys.stderr)
            print(mock_stderr.getvalue(), file=sys.stderr)
            print("STDOUT WAS:", file=sys.stderr)
            print(mock_stdout.getvalue(), file=sys.stderr)
            raise
    
    def assertEnergyPrinted(self, stdout, *energies):
        stdout_str=stdout.getvalue()
        # Energy printed to stdout
        self.assertIn("total_energy", stdout_str)
        # The default energy consists of ROG, SLD and AME energies.
        for energy in energies:
            self.assertIn(energy, stdout_str)
    
    def getStatsFor(self, stdout, name, withA=False):
        stdout.seek(0)
        header = None
        data=[]
        for line in stdout:
            if header is None:
                if not line.strip().startswith("#"):
                    header=line.split("\t")                
            else:
                fields = line.strip().split("\t") #The last line!

                val = fields[header.index(name)]
                if withA:
                    val=val.split()[0]
                data.append(float(val))
        return data
    def getSavedFile(self, mock, filename):
        for i_rev, call in enumerate(reversed(mock.mock_calls)):
            if call[1]:
                if call[1][0] == filename:
                    break
        else:
            return None
        for i in range(-i_rev, len(mock.mock_calls)):
            call = mock.mock_calls[i]
            if call[0]=="().write":
                cg = ftmc.CoarseGrainRNA()
                cg.from_cg_string(call[1][0])
                return cg
        else:
            return None
    #def printFilenames(self, mock):
    #    for call in mock.mock_calls:
    #        if not call[0]:
    #            print(call[1])
        
    def allSavedFiles(self, mock, pattern):
        cgs=[]
        in_file = False
        for call in mock.mock_calls:
            if call[1]:
                try:
                    if pattern.match(call[1][0]):
                        in_file=True
                except: pass
            if in_file and call[0]=="().write":
                cg = ftmc.CoarseGrainRNA()
                cg.from_cg_string(call[1][0])
                cgs.append(cg)
                in_file=False
        return cgs
    def assertStatsOnlyDifferFor(self, cgs, element):
        prev_stats={}
        for elem in cgs[0].defines:
            if elem[0] in "mi":
                prev_stats[elem]=cgs[0].get_bulge_angle_stats(elem)
        for i, cg in enumerate(cgs[1:]):
            for elem in cg.defines:
                if elem!=element and elem[0] in "mi":
                    self.assertEqual(cg.get_bulge_angle_stats(elem), prev_stats[elem],
                                           msg="Changesd stats for bulge {} in step {}".format(elem, i+1))
    def countDifferentStatsFor(self, cgs, element):
        stats=set()
        for cg in cgs:
            if element[0] in "mi":
                stat = cg.get_bulge_angle_stats(element) #Hash and eq should be implemented properly
            elif element[0] =="s":
                stat = str(cg.get_stem_stats(element))
            else:
                stat= str(cg.get_loop_stat(element))
            stats.add(stat)
        return len(stats)

    def assertManyDifferentStatsFor(self, cgs, element, n):
        self.assertGreater(self.countDifferentStatsFor(cgs, element), n)

class ErnwinTestBase(unittest.TestCase, ErnwinTestsMixin):
    def setUp(self):
        self.parser = ernwin.get_parser()
        self.longMessage = True
        # We need to reset the stats in both setup and tearDown, 
	# because tests from other files are not guaranteed to do the same.
        ftms.ConstructionStats.angle_stats = None
        ftms.ConstructionStats.stem_stats = None
        ftms.ConstructionStats.loop_stats = None
        ftms.ConstructionStats.fiveprime_stats = None
        ftms.ConstructionStats.threeprime_stats = None
        ftms.ConstructionStats.conf_stats = None
    def tearDown(self):
        ftms.ConstructionStats.angle_stats = None
        ftms.ConstructionStats.stem_stats = None
        ftms.ConstructionStats.loop_stats = None
        ftms.ConstructionStats.fiveprime_stats = None
        ftms.ConstructionStats.threeprime_stats = None
        ftms.ConstructionStats.conf_stats = None

class TestCommandLineUtilGeneralBehaviour(ErnwinTestBase):
    def test_start_from_scratch(self):
        open_main, open_stats, stdout, stderr = self.runErnwin(
                              "python ernwin_new.py test/fess/data/1GID_A-structure1.coord -i 1 "
                              "--step-save 1 --start-from-scratch --seed 1")
        orig_cg = ftmc.CoarseGrainRNA("test/fess/data/1GID_A-structure1.coord")
        new_cg = self.getSavedFile(open_stats, "1GID_A/step000001.coord")
        self.assertGreater(ftmsim.cg_rmsd(orig_cg, new_cg), 5)

    def test_start_not_from_scratch(self):
        open_main, open_stats, stdout, stderr = self.runErnwin(
                              "python ernwin_new.py test/fess/data/1GID_A-structure1.coord -i 1 "
                              "--step-save 1 --seed 1")
        orig_cg = ftmc.CoarseGrainRNA("test/fess/data/1GID_A-structure1.coord")
        new_cg = self.getSavedFile(open_stats, "1GID_A/step000001.coord")
        self.assertLess(ftmsim.cg_rmsd(orig_cg, new_cg), 5)

    def test_exhaustive(self):
        open_main, open_stats, stdout, stderr = self.runErnwin(
                              "python ernwin_new.py test/fess/data/1GID_A-structure1.coord -i 50 "
                              "--step-save 1 --exhaustive i0 --seed 1")
        cgs = self.allSavedFiles(open_stats, re.compile("1GID_A\/step.*\.coord"))
        for cg in cgs:
            self.assertStatsOnlyDifferFor(cgs, "i0")
        self.assertManyDifferentStatsFor(cgs, "i0", 40)

    def test_new_ml(self):
        #We use only ROG energy because it is faster than the A-Minor energy
        
        #Old ML sampling
        open_main, open_stats, stdout, stderr = self.runErnwin(
                              "python ernwin_new.py test/fess/data/4way.cg -i 2000 "
                              "--step-save 1 --seed 1 --energy ROG")
        normal_cgs = self.allSavedFiles(open_stats, re.compile(".*\/step.*\.coord"))

        # New ML sampling
        open_main, open_stats, stdout, stderr = self.runErnwin(
                              "python ernwin_new.py test/fess/data/4way.cg -i 2000 "
                              "--step-save 1 --new-ml --seed 1 --energy ROG")
        newml_cgs = self.allSavedFiles(open_stats, re.compile(".*\/step.*\.coord"))

        #Count how many different stats per multiloop segment
        normal_counts = [self.countDifferentStatsFor(normal_cgs, m) for m in ["m0", "m1", "m2", "m3"]]
        new_counts = [self.countDifferentStatsFor(newml_cgs, m) for m in ["m0", "m1", "m2", "m3"] ]

        if max(normal_counts)/min(normal_counts)<5:
            raise WrongAssumptionAboutTestInput("With old multiloop we expect to find a different "
                                                "number of stats per multiloop segment.")                                                
        # Note that the multiloop in 4way.cg has a symmetrical bulge graph (except for a missing 
        # hairpin at the 3'/5' stem). For this reason, perfect sampling should change each 
        # multiloop segment equally often.
        #print(normal_counts)
        #print(new_counts)
        self.assertLess(max(new_counts)/min(new_counts), 1.5, msg="The ratio "
                  "max(number of stats per element)/min(number of stats per element) is too high."
                  " We sample different multiloop segments not equally!")

class TestCommandLineUtilOutputOptions(ErnwinTestBase):
    def test_dist(self):
        open_main, open_stats, stdout, stderr = self.runErnwin(
                "python ernwin_new.py test/fess/data/1GID_A-structure1.coord -i 100 "
                "--dist 1,20:2,21 --seed 1")
        print(stdout.getvalue())
        d1 = np.array(self.getStatsFor(stdout, "Distance_1-20", True))
        d2 = np.array(self.getStatsFor(stdout, "Distance_2-21", True))
        np.testing.assert_almost_equal(d1,d2, decimal=-1)

class TestCommandLineUtilStats(ErnwinTestBase):
    def test_empty_stats_file(self):
        with self.assertRaises(LookupError):
            open_main, open_stats, stdout, stderr = self.runErnwin(
                "python ernwin_new.py test/fess/data/1GID_A-structure1.coord -i 100 "
                "--stats-file test/fess/data/empty.stats --seed 1")
    def test_statsfile_single_stat_per_elem(self):
        open_main, open_stats, stdout, stderr = self.runErnwin(
            "python ernwin_new.py test/fess/data/4way.cg -i 1 --step-save 1 "
            "--stats-file test/fess/data/statsFor4way.stats --seed 1")
        cgs = self.allSavedFiles(open_stats, re.compile(".*\/step.*\.coord"))
        new_counts = [self.countDifferentStatsFor(cgs, d) for d in cgs[0].defines ]
        self.assertLess(max(new_counts),3) #Two different stats are possible: The original and the one from the file.
    @unittest.skip("Clustered angle stats do not yet work with the new stats container")
    def test_clustered_angle_stats(self):
        open_main, open_stats, stdout, stderr = self.runErnwin(
            "python ernwin_new.py test/fess/data/4way.cg -i 500 --step-save 1 "
            "--clustered-angle-stats fess/stats/clustered_angular_0.3.stats --seed 1")
        anglestats = ftms.ClusteredAngleStats("fess/stats/clustered_angular_0.3.stats")
        print(anglestats, type(anglestats))
        cgs = self.allSavedFiles(open_stats, re.compile(".*\/step.*\.coord"))
        cgs[0].traverse_graph()
        if "m2" not in cgs[0].mst:
            raise WrongAssumptionAboutTestInput("m2 expected to be in mst")
        clusters = []
        oldstat = None
        for cg in cgs:
            stat = cg.get_bulge_angle_stats("m2")
            if stat != oldstat: #Only when this element was resampled
                clusters.append(anglestats.cluster_of(stat[0]))
                oldstat = stat
        if len(clusters)<4:
            raise WrongAssumptionAboutTestInput("Too high reject rate for m2. Expecting to see at least 4 different stats.")
        self.assertEqual(len(clusters), len(set(clusters)))
    @unittest.skip("Jar3d does not work yet with new stats container")
    def test_jar3d(self):
        cg = ftmc.CoarseGrainRNA("test/fess/data/1GID_A-structure1.coord")
        open_main, open_stats, stdout, stderr = self.runErnwin(
                "python ernwin_new.py test/fess/data/1GID_A-structure1.coord -i 100 "
                "--jar3d --step-save 1")
        cgs = self.allSavedFiles(open_stats, re.compile(".*\/step.*\.coord"))
        self.assertEqual(self.countDifferentStatsFor(cgs, "i4"), 2, msg="2 stats: the original and only 1 jar3d hit")

def patchSeed():
    original_seed = random.seed
    def f(*args):
        print("SEEDING RANDOM WITH {}".format(args), file=sys.stderr)
        return original_seed(*args)
    return f


class TestCommandLineUtilEnergyDefault(ErnwinTestBase):
    def test_default(self):
        if True: #with patch('random.seed', side_effect=patchSeed(), autospec=True):
            open_main, open_stats, stdout, stderr = self.runErnwin(
                              "python ernwin_new.py test/fess/data/1GID_A-structure1.coord -i 10 --seed 1"
                                                              )

        # Ernwin with default configuration creates at least the following files
        # (The folder name depends on the Name in the coord/ fasta-file of the initial structure):
        open_main.assert_any_call('1GID_A/out.log', 'w')
        open_stats.assert_any_call('1GID_A/best0.coord', 'w')
        open_stats.assert_any_call('1GID_A/best_rmsd0.coord', 'w')
        open_stats.assert_any_call('1GID_A/best1.coord', 'w')
        open_stats.assert_any_call('1GID_A/best_rmsd1.coord', 'w')
        open_stats.assert_any_call('1GID_A/best2.coord', 'w')
        open_stats.assert_any_call('1GID_A/best_rmsd2.coord', 'w')
        # The -i option defines the number of iterations. 
        # The out-file contains 1 line for the header and some comments plus one line per iteration
        self.assertOneWritePerIteration(open_main().write, 10)

class TestCommandLineUtilEnergyOption(ErnwinTestBase):
    def test_energy_evaluation_DEF(self):
        open_main, open_stats, stdout, stderr = self.runErnwin(
                      "python ernwin_new.py test/fess/data/1GID_A-structure1.coord --eval-energy"
                                                              )
        # No files written in energy evaluation mode (at least asserting this for some modules)
        self.assertEqual(open_main.call_count, 0)
        self.assertEqual(open_stats.call_count, 0)
        # Energy printed to stdout
        
        self.assertEnergyPrinted(stdout, "ROG", "AME(0)", "AME(1)", "SLD")

    def test_energy_evaluation_NDR(self):
        open_main, open_stats, stdout, stderr = self.runErnwin(
                                    "python ernwin_new.py test/fess/data/1GID_A-structure1.coord "
                                    "--eval-energy --energy NDR20.2"
                                                              )
        # No files written in energy evaluation mode (at least asserting this for some modules)
        self.assertEqual(open_main.call_count, 0)
        self.assertEqual(open_stats.call_count, 0)
        # Energy printed to stdout
        self.assertIn("total_energy", stdout.getvalue(), msg=stdout.getvalue())
        self.assertEnergyPrinted(stdout, "NDR20.2")

    def test_NDR_energy_can_change_rog(self):        
        #Originally, the ROG is small:
        origROG = ftmc.CoarseGrainRNA("test/fess/data/1GID_A-structure1.coord").radius_of_gyration()
        if origROG > 30.0:
            raise WrongAssumptionAboutTestInput("1GID_A-structure1.coord should have a smaller "
                                                "initial ROG. Found {}".format(origROG))
        # The NDR60 energy will try to generate a radius of gyration of ~42
        open_main, open_stats, stdout, stderr = self.runErnwin(
                                     "python ernwin_new.py test/fess/data/1GID_A-structure1.coord "
                                     "-i 500 --energy NDR60 --seed 1")
        rog = self.getStatsFor(stdout, "ROG", True)[-1]
        self.assertGreater(float(rog), 40.0)

    def test_CHE_energy_reduces_rmsd(self):
        open_main, open_stats, stdout, stderr = self.runErnwin(
                                     "python ernwin_new.py test/fess/data/1GID_A-structure1.coord "
                                     "-i 500 --energy CHE --seed 1 --start-from-scratch")
        stdout.seek(0)
        firstRMSD=None
        header = None
        for line in stdout:                
            if line.strip().startswith("#"): continue
            if header is None:
                header=line.split("\t")
            elif firstRMSD is None:
                fields = line.strip().split("\t") #The last line!
                try:
                    firstRMSD = fields[header.index("RMSD")].split()[0]
                except ValueError:
                    print("HEADER:", header)
                    raise
        fields = line.strip().split("\t") #The last line!
        lastRMSD = float(fields[header.index("RMSD")].split()[0])
        self.assertGreater(firstRMSD, lastRMSD-100)

    def test_CLA_energy_clamps_together(self):
        #Make sure the original distances are not desired
        orig_cg = ftmc.CoarseGrainRNA("test/fess/data/1GID_A-structure2.coord")
        h0h2_orig = orig_cg.element_physical_distance("h0","h2")
        h1h2_orig = orig_cg.element_physical_distance("h0","h1")
        print(h0h2_orig, h1h2_orig)
        if h1h2_orig>h0h2_orig-10:
            raise WrongAssumptionAboutTestInput("In the initial cg model, h2 should be "
                    "significantly closer to h0 than to h1. Found h0-h2={}, "
                    "h1-h2={}".format(h0h2_orig, h1h2_orig))

        open_main, open_stats, stdout, stderr = self.runErnwin(
                "python ernwin_new.py test/fess/data/1GID_A-structure1.coord "
                "-i 500 --energy CLA --clamp h0,h2  --seed 1 --step-save 500")
        new_cg = self.getSavedFile(open_stats, "1GID_A/step000500.coord")
        h0h2 = new_cg.element_physical_distance("h0","h2")
        h1h2 = new_cg.element_physical_distance("h0","h1")
        print(h0h2, h1h2)
        self.assertLess(h0h2, h1h2)
        self.assertLess(h0h2, 30)
        self.assertLess(h0h2, h0h2_orig)
    @unittest.skip("We need to fix the Projection based stuff")
    def test_PRO_energy(self):
        # We use distances from a projection of 1GID_A_structure2.coord, 
        # proj.dir (-0.147,-0.311,-0.876)
        open_main, open_stats, stdout, stderr = self.runErnwin(
              "python ernwin_new.py test/fess/data/1GID_A-structure1.coord "
              "-i 500 --energy 50PRO --seed 1 --step-save 500 --projected-dist h0,h1,31.64:"
              "h0,h2,13.56:h0,s3,58.44:h1,h2,30.30:h1,s3,54.27:h2,s3,44.95")
        # Although we start from structure1, our result should stronger resemble structure2,
        # Because we used distances from structure2
        orig_cg = ftmc.CoarseGrainRNA("test/fess/data/1GID_A-structure1.coord")
        target_cg = ftmc.CoarseGrainRNA("test/fess/data/1GID_A-structure2.coord")
        new_cg = self.getSavedFile(open_stats, "1GID_A/step000500.coord")
        self.assertLess(ftmsim.cg_rmsd(new_cg, target_cg), ftmsim.cg_rmsd(new_cg, orig_cg))

    # This test is very slow!
    @unittest.skip("This test takes too long")
    def test_HDE_energy(self):
        open_main, open_stats, stdout, stderr = self.runErnwin(
              "python ernwin_new.py test/fess/data/1GID_A-structure1.coord "
              "-i 200 --energy HDE --seed 1 --step-save 200 "
              "--ref-img test/fess/data/1GID_A-structure2.coord.dpi15.width120.4.png --scale 120")
        # Although we start from structure1, our result should stronger resemble structure2,
        # Because we used distances from structure2
        orig_cg = ftmc.CoarseGrainRNA("test/fess/data/1GID_A-structure1.coord")
        target_cg = ftmc.CoarseGrainRNA("test/fess/data/1GID_A-structure2.coord")
        new_cg = self.getSavedFile(open_stats, "1GID_A/step000200.coord")
        self.assertLess(ftmsim.cg_rmsd(new_cg, target_cg), ftmsim.cg_rmsd(new_cg, orig_cg))

    # This test is very slow!
    @unittest.skip("This test takes too long")
    def test_FPP_energy(self):
        open_main, open_stats, stdout, stderr = self.runErnwin(
              "python ernwin_new.py test/fess/data/1GID_A-structure1.coord "
              "-i 200 --energy FPP --seed 1 --step-save 200 --scale 120 "
              "--fpp-landmarks 21,7,12:1,10,7:137,8,4:50,7,2 --ref-img test/fess/data/1GID_A-structure2.forFPP.png")

        # Although we start from structure1, our result should stronger resemble structure2,
        # Because we used distances from structure2
        orig_cg = ftmc.CoarseGrainRNA("test/fess/data/1GID_A-structure1.coord")
        target_cg = ftmc.CoarseGrainRNA("test/fess/data/1GID_A-structure2.coord")
        new_cg = self.getSavedFile(open_stats, "1GID_A/step000200.coord")
        self.assertLess(ftmsim.cg_rmsd(new_cg, target_cg), ftmsim.cg_rmsd(new_cg, orig_cg))

class TestCombinedOptions(ErnwinTestBase):
    def test_clamp_with_dist(self):
        cg = ftmc.CoarseGrainRNA("test/fess/data/1GID_A-structure1.coord")
        open_main, open_stats, stdout, stderr = self.runErnwin(
                "python ernwin_new.py test/fess/data/1GID_A-structure1.coord -i 200 "
                "--dist 65,136 --seed 1 --energy ROG,SLD,AME,CLA --clamp 65,136")#Clamping h1 and h2

        dists = np.array(self.getStatsFor(stdout, "Distance_65-136", True))
        print(dists)
        
        self.assertLess(np.mean(dists[100:]), np.mean(dists[:100])) #Decreasing dist
        self.assertLess(np.mean(dists[150:]), 25) # Minimal clamp energy is at 15A
        #self.assertLess(np.mean(dists[180:]), 18) # Minimal clamp energy is at 15A
    @unittest.skip("Jar3d does not work yet with new stats container")
    def test_jar3d_with_exhaustive(self):
        cg = ftmc.CoarseGrainRNA("test/fess/data/1GID_A-structure1.coord")
        open_main, open_stats, stdout, stderr = self.runErnwin(
                "python ernwin_new.py test/fess/data/1GID_A-structure1.coord -i 10 "
                "--jar3d --exhaustive i4")
        print (stdout.getvalue())
        rogs = np.array(self.getStatsFor(stdout, "ROG", True))
        self.assertEqual(max(rogs[1:]), min(rogs[1:]))
