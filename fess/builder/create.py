"""
A module for functions used during structure creation,
either via sampling or via building.

This module contains functions used both by builder.py and move.py
"""
from __future__ import absolute_import, division, print_function
from builtins import (ascii, bytes, chr, dict, filter, hex, input, #pip install future
                      int, map, next, oct, open, pow, range, round,
                      super, zip) # Do not import str here (incompatible with pyffx)!
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)
__metaclass__=type

import math
import random
import itertools as it
import logging



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
    try:
        import pyffx
    except:
        warnings.warn("Could not import pyffx. Sampling stat combinations in a "
                      "slower/ memory inefficient way.")
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
    else:
        log.info("Sampling stat combinations via format preserving encryption.")
        for value in _random_iteration(choices):
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

def index_to_sample(index, choices):
    sampled = {}
    for elem, stats in choices.items():
        c_nr = index%len(stats)
        sampled[elem]=stats[c_nr]
        index = int(index/len(stats))
    return sampled


def _random_iteration(choices):
    """
    Use format preserving encryption (via pyffx) to iterate through all
    combinations of stats in a pseudorandom order.

    Note: pyffx allows setting of the number of rounds and is (thus) faster than libffx
    """
    import pyffx
    product_size = 1
    for stats in choices.values():
        product_size *= len(stats)
    # The pyffx does not allow setting a max-value but instead
    # a (maximal) number of digits. For this reasons we convert to binary
    # and reject everything that Format preserving encryption maps to a too high
    # value. In the case of decimal, we would expect 9/10 of all numbers to
    # be rejected, while with binary only half of all numbers are rejected.
    string_length = (int(math.log(product_size, 2))+1)
    log.debug("product_size: %s, string_length %s", product_size, string_length)
    fpe = pyffx.String(pyffx.FFX(str(random.random()),rounds=5), alphabet="01", length=string_length)
    format_string = "{:0"+str(string_length)+"b}"
    for i in range(2**string_length):
        log.debug("encrypting (%s) %r", format_string, format_string.format(i))
        perm_number_binary = fpe.encrypt(format_string.format(i))
        perm_number = int(perm_number_binary, 2)
        log.debug("FPE mapping %d -> %d (p_size=%d)", i, perm_number, product_size)
        if perm_number<product_size:
            yield index_to_sample(perm_number, choices)
        else:
            log.debug("Discarding %d", perm_number)
