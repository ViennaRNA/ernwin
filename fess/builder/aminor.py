from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input, #pip install future
                      map, next, oct, open, pow, range, round,
                      str, super, zip)
__metaclass__=object

try:
    from collections.abc import Set
except ImportError:
    from collections import Set

import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.vector as ftuv
import forgi.utilities.debug as fud
import forgi.graph.bulge_graph as fgb
from collections import namedtuple
import warnings
import logging
import numpy as np
import scipy.stats
import os.path

from logging_exceptions import log_to_exception

log=logging.getLogger(__name__)

try:
  profile  #The @profile decorator from line_profiler (kernprof)
except:
  def profile(x):
    return x


def total_prob(individual_probs):
    """
    Return the total probability for the loop to form at least one A-Minor interaction.

    This function tries for interaction of the loop with
    all stems (except adjacent ones) and returns p1+(1-p1)*p2+...
    where p1, p2, ... are the individual probabilities.
    """
    log.debug("Entering 'total_prob'")
    total_prob = 0
    for p in individual_probs:
        log.debug("Adding p = %f to total_prob = %f", p, total_prob)
        total_prob += (1-total_prob)*p
    log.debug("total_prob: Returning: %s", total_prob)
    if total_prob>1:
        log.error("Probability >1 for %s %s with domain %s", cg.name, loop, domain)
        assert False
    return total_prob

def iter_probs(loop, cg, prob_fun, cutoff_dist, domain = None):
    """
    Iterate over all stems and yield the probability for an interaction with loop.

    .. warning ::

        Do not rely on len(list(_iter_probs))!
        This function yields a single zero as a last element,
        which does not correspond to any stem.
        This is required to make sure at least one value is yielded in cases
        where no stem is close enough for interactions.

    For params: See the documentation of total_prob
    """
    log.debug("Entering '_iter_probs'")
    if domain is not None:
        stems = (s for s in domain if s[0]=="s")
    else:
        stems = cg.stem_iterator()
    for s in stems:
        if s in cg.edges[loop]:
            continue
        if not ftuv.elements_closer_than(cg.coords[loop][0],
                                 cg.coords[loop][1],
                                 cg.coords[s][0],
                                 cg.coords[s][1],
                                 cutoff_dist):
            continue
        point = get_relative_orientation(cg, loop, s)
        p, = prob_fun(point)
        if p>1:
            warnings.warn("Probability at %s is %f>1 for %s %s with domain %s" %(point, p, cg.name, loop, domain))
        yield(p)
    yield 0 #Always yield at least one number, so max(_iter_probs(...)) does not raise an error
