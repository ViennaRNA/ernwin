#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input, #pip install future
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)
__metaclass__=object

class Mover:
    def __init__(self, stat_source):
        self.stat_source = stat_source
        #: A tuple (elemenmt_name, stat)
        self._prev_stat = None
    def move(self, sm):
        """
        Propose a single Monte Carlo move.
        
        :param sm: The spatial model
        """
        possible_elements = sm.bg.get_mst()
        d = random.choice(list(possible_elements))
        possible_stats = self.stat_source.get_possible_stats(sm.bg, d)
        new_stat = random.choice(possible_stats)
        self._prev_stat = (d, sm.elem_defs[d])
        sm.elem_defs[d]=new_stat
        sm.new_traverse_and_build(start = d, include_start = True)
        movestring = "{}:{}->{};".format(d, self._prev_stat[1].pdb_name, new_stat.pdb_name)
        return movestring
        
    def revert(self, sm):
        """
        Revert the last Monte Carlo move performed by this mover.
        """
        pass
