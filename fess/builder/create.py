"""
A module for functions used during structure creation,
either via sampling or via building.

This module contains functions used both by builder.py and move.py
"""
from __future__ import absolute_import, division, print_function
from builtins import (ascii, bytes, chr, dict, filter, hex, input, #pip install future
                      int, map, next, oct, open, pow, range, round,
                      super, zip) # Do not import str here (incompatible with pyffx)!
__metaclass__=type

import math
import random
import itertools as it
import logging

has_warned = set()

log = logging.getLogger(__name__)

def stat_combinations(cg, elems, stat_source, unique = False):
    """
    Yield combinations of stats for the elements.

    If the total number of combinations is large, this samples combinations
    (with the possibility of duplicates). If the total number of
    combinations is small, this iterates through all combinations in a
    pseudorandom order (no duplicates).
    If the number of combinations is large but not huge,
    """
    choices = { elem : list(stat_source.iterate_stats_for(cg, elem))
                for elem in elems }
    for value in _sample_unique_stat_combinations(choices):
        yield value

def index_to_sample(index, choices):
    sampled = {}
    for elem, stats in choices.items():
        c_nr = index%len(stats)
        sampled[elem]=stats[c_nr]
        index = int(index/len(stats))
    return sampled

def _sample_unique_stat_combinations(choices):
    from bitarray import bitarray
    product_size = 1
    for stats in choices.values():
        product_size *= len(stats)
    if product_size>10**9:
        if tuple(choices.keys()) not in has_warned:
            log.warning("Too many stats-combination for unique sampling of %s: %s."
                        "Sampling with replacement.", list(choices.keys()), product_size)
            has_warned.add(tuple(choices.keys()))
        while True:
            i = random.randrange(product_size)
            yield index_to_sample(i, choices)
    else:
        found = bitarray(product_size)
        found.setall(False)
        while not found.all():
            i = random.randrange(product_size)
            if found[i]:
                continue
            found[i]=True
            yield index_to_sample(i, choices)
