from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)
from future.utils import viewkeys
import unittest
import sys
import StringIO
import fess.builder.stat_container as fbstat
import forgi.threedee.model.stats as ftmstats

from mock import mock_open, patch

class ParseFileTests(unittest.TestCase):
    def test_parse_empty_line(self):
        filecontent = StringIO.StringIO("\nstem test:s_0 5 10.388 2.43294047108")
        stats = fbstat.parse_stats_file(filecontent)
        self.assertEqual(len(stats["stem"]), 1)
        self.assertEqual(len(stats["angle"]), 0)
        self.assertEqual(stats["stem"][5][0].bp_length, 5)
        
    def test_parse_line_with_comments(self):
        filecontent = StringIO.StringIO("\nstem test:s_0 5 10.388 2.43294047108\n"
                                        " # Das ist ein comment\n"
                                        "angle test:i_0 5 2 2.203691 2.099941 0.586450 17.134279 1.191397 1.274896 1\n"
                                        "angle test:i_1 4 2 2.203691 2.099941 0.586450 17.134279 1.191397 1.274896 1 #Another # comment")
        stats = fbstat.parse_stats_file(filecontent)
        self.assertEqual(len(stats["stem"]), 1)
        self.assertEqual(len(stats["angle"]), 2)
        self.assertEqual(stats["angle"][(5,2,1)][0].pdb_name, "test:i_0")
        self.assertEqual(stats["angle"][(4,2,1)][0].pdb_name, "test:i_1")

    

class ReadFileTests(unittest.TestCase):
    def setUp(self):
        pass
    def test_read_stats_file(self):
        stats = fbstat.read_stats_file("test/fess/data/test1.stats")
        self.assertEqual(len(stats["stem"]), 1)
        self.assertEqual(len(stats["angle"]), 3)
        self.assertEqual(len(stats["loop"]), 1)
        self.assertEqual(len(stats["3prime"]), 1)
        self.assertEqual(len(stats["5prime"]), 1)
        
        self.assertEqual(stats["stem"][5], 
                         [ftmstats.StemStat("stem test:s_0 5 10.388 2.43294047108")])
        a_stat = ftmstats.AngleStat()
        a_stat.parse_line("angle test:i_0 5 2 2.203691 2.099941 0.586450 17.134279 1.191397 1.274896 1")
        self.assertEqual(stats["angle"][(5, 2, 1)], 
                         [a_stat])
        a_stat.parse_line("angle test:m_0 4 1000 0.985166 -1.185545 -2.000463 13.701389 0.982669 0.267821 -4")
        self.assertEqual(stats["angle"][(4, 1000, -4)], 
                         [a_stat])
        self.assertEqual(stats["loop"][6],
                         [ftmstats.LoopStat("loop test:h_0 6  15.2401560955 0.269833051418 0.731484795668")])
        self.assertEqual(stats["3prime"][4],
                         [ftmstats.LoopStat("3prime test:t_0 4 20.4034805163 1.47912394946 -0.0715301558972")])
        self.assertEqual(stats["5prime"][4],
                         [ftmstats.LoopStat("5prime test:f_0 4 20.4034805163 1.47912394946 -0.0715301558972")])

class StatStorageTest(unittest.TestCase):      
    def test_stat_files_are_loaded_lazily(self):        
        stats_open = mock_open()
        with patch('fess.builder.stats_container', stats_open, create=True):
            st = fbstat.StatStorage("test/fess/data/test1.stats", ["test/fess/data/fallback1.stats", "test/fess/data/fallback2.stats"])
        self.assertEqual(len(stats_open.mock_calls), 0) 

    def test__iter_stat_sources(self):
        st = fbstat.StatStorage("test/fess/data/test1.stats", ["test/fess/data/fallback1.stats", "test/fess/data/fallback2.stats"])
        source_iter = st._iter_stat_sources()
        first = next(source_iter)
        self.assertEqual(viewkeys(first["stem"]), {5})
        self.assertEqual(len(first["stem"][5]), 1)
        second = next(source_iter)
        self.assertEqual(viewkeys(second["stem"]), {6})
        self.assertEqual(len(second["stem"][6]), 1)
        third = next(source_iter)
        self.assertEqual(viewkeys(third["stem"]), {5, 10})
        self.assertEqual(len(third["stem"][5]), 2)

    def test__possible_stats_without_fallback(self):
        st = fbstat.StatStorage("test/fess/data/test1.stats")
        try:
            stat = st._possible_stats("stem", 5)[0]
        except:
            print("SOURCES", st._sources)
            raise
        self.assertEqual(stat.pdb_name, "test:s_0")
        with self.assertRaises(LookupError):
            st._possible_stats("stem", 12) #No stem with length 12
        self.assertEqual(st._possible_stats("angle", (5,2,1))[0].pdb_name, "test:i_0")
        with self.assertRaises(LookupError):
            st._possible_stats("angle", (4,5,1))
        self.assertEqual(st._possible_stats("angle", (4,1000,-4))[0].pdb_name, "test:m_0")
        with self.assertRaises(LookupError):
            st._possible_stats("angle", (4,1000,3))
        self.assertEqual(st._possible_stats("loop", 6)[0].pdb_name, "test:h_0")
        with self.assertRaises(LookupError):
            st._possible_stats("loop", 4)
        self.assertEqual(st._possible_stats("3prime", 4)[0].pdb_name, "test:t_0")
        with self.assertRaises(LookupError):
            st._possible_stats("3prime", 6)
        self.assertEqual(st._possible_stats("5prime", 4)[0].pdb_name, "test:f_0")
        with self.assertRaises(LookupError):
            st._possible_stats("5prime", 6)

    def test_pick_stat_with_fallback(self):
        st = fbstat.StatStorage("test/fess/data/test1.stats", ["test/fess/data/fallback1.stats", "test/fess/data/fallback2.stats"])
        #If the first file has enough stats, do not fall back
        poss_stats = st._possible_stats("stem", 5, 1)
        self.assertEqual(poss_stats[0].pdb_name, "test:s_0")
        self.assertEqual(len(poss_stats), 1)
        #No hits in first file, fall back to second
        poss_stats = st._possible_stats("stem", 6, 1)
        self.assertEqual(poss_stats[0].pdb_name, "fallback1:s_0")
        self.assertEqual(len(poss_stats), 1)
        #No hits in second file, fall back to third
        poss_stats = st._possible_stats("stem",10, 1)
        self.assertEqual(poss_stats[0].pdb_name, "fallback2:s_1")
        self.assertEqual(len(poss_stats), 1)
        #We want at least two hits. Collect them from 1st and 3rd file
        poss_stats = st._possible_stats("stem",5, 2)
        self.assertEqual(len(poss_stats), 3) #Even if we requested a minimum of 2 stats, we end up with 3.















