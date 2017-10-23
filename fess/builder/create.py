"""
A module for functions used during structure creation,
either via sampling or via building.

This module contains functions used both by builder.py and move.py
"""
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input, #pip install future
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)
__metaclass__=type


import random
import itertools as it


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
    target_len = 1
    for elem in elems:
        target_len*=len(choices[elem])
    if target_len < 5*10**6: #list(range(5*10^6)) should fit in 500MB Memory
        stat_gen = _iter_stat_combinations
    else:
        stat_gen = _sample_stat_combinations
    seen = set()
    for value in stat_gen( choices):
        if unique:
            if value in seen:
                continue
            seen.add(value)
        yield value

def _sample_stat_combinations(choices):
    elems = list(choices.keys())
    max_i = max( len(stats) for stats in choices.values())
    while True:
        sampled = {}
        for elem in elems:
            try:
                sampled[elem] = choices[elem][random.randrange(max_i)]
            except IndexError:
                # If one index is out of range, discard the whole sample and try again
                # This is necessary so all combinations have equal probability
                break
        else: # no break
            yield sampled

def _iter_stat_combinations(choices):
    """
    Iterate over combinations of stats in a pseudorandom order.
    This only works for small numbers of elems/ stats.

    :param cg: A CoarseGrainRNA object
    :param elem: A list of element names
    :param stat_source: A instance of fess.builder.stat_container.StatStorage
    :yields: A dictionary {elem:stat}
    """
    elems = list(choices.keys())
    indices = list(it.product(*(range(len(choices[elem])) for elem in elems)))
    random.shuffle(indices)
    for v in choices.values():
        random.shuffle(v)
    for index_combi in indices:
        sampled = {}
        for i, elem in enumerate(elems):
            sampled[elem] = choices[elem][index_combi[i]]
        yield sampled
