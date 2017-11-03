import unittest
import logging


import fess.builder.create as fbc

log = logging.getLogger(__name__)

class TestFFX(unittest.TestCase):
    def setUp(self):
        self.choices = {"lc":"abcde", "uc": "ABC", "num": "12", "sym": "!@#$%^&"}
    def test_runs(self):
        a=next(fbc._random_iteration(self.choices))
        assert "lc" in a

    def test_all_values_produces_once(self):
        all_sampled = set()
        for choice in fbc._random_iteration(self.choices):
            choice = "{lc}{uc}{num}{sym}".format(**choice)
            self.assertNotIn(choice,  all_sampled)
            all_sampled.add(choice)
        log.error(all_sampled)
        self.assertEqual(len(all_sampled), 210)
        self.assertIn("aB2%", all_sampled)

class TestBitarray(unittest.TestCase):
    def setUp(self):
        self.choices = {"lc":"abcde", "uc": "ABC", "num": "12", "sym": "!@#$%^&"}
    def test_runs(self):
        a=next(fbc._sample_stat_combinations(self.choices))
        assert "lc" in a

    def test_all_values_produces_once(self):
        all_sampled = set()
        for choice in fbc._sample_unique_stat_combinations(self.choices):
            choice = "{lc}{uc}{num}{sym}".format(**choice)
            self.assertNotIn(choice,  all_sampled)
            all_sampled.add(choice)
        log.error(all_sampled)
        self.assertEqual(len(all_sampled), 210)
        self.assertIn("aB2%", all_sampled)
