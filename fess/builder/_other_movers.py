from __future__ import absolute_import
from __future__ import print_function
from builtins import zip

import random
import itertools as it
import copy
import logging

from .move import Mover
from six.moves import range

log=logging.getLogger(__name__)


"""
Implementations of other Mover subclasses, which do not consistently perform well
in benchmarks, but might be usefull in certain special applications.
"""

class MSTchangingMover(Mover):
    def __init__(self, stat_source, **kwargs):
        super(MSTchangingMover, self).__init__(stat_source)
        self._prev_mst = None
    def move(self, sm):
        self._prev_stats = {}
        self._prev_mst = None
        # Get an element. If it is a ml-segment, break it.
        elem = self._get_elem(sm)
        if elem[0]=="m":
            loop = sm.bg.shortest_mlonly_multiloop(elem)
            if "open" in sm.bg.describe_multiloop(loop):
                return super(MSTchangingMover, self).move(sm)
            else:
                defined_loop = set(loop)&sm.bg.mst
                self._prev_mst = copy.copy(sm.bg.mst)
                for node in defined_loop:
                    self._prev_stats[node]= sm.elem_defs[node]
                old_elem = elem
                # Now change the MST
                elem = sm.set_multiloop_break_segment(elem)
                if elem is None:
                    # The element was already broken => call super
                    return super(MSTchangingMover, self).move(sm)
                assert elem.startswith("m")
                new_stat = self.stat_source.sample_for(sm.bg, elem)
                sm.elem_defs[elem]=new_stat
                sm.new_traverse_and_build(start = elem, include_start = True)
                return "BREAK{};{}->{}".format(old_elem, elem, new_stat.pdb_name)
        else:
            return super(MSTchangingMover, self).move(sm)

    def revert(self, sm):
        if self._prev_mst is not None:
            # Reset MST
            log.info("Reverting MST")
            sm.change_mst(self._prev_mst, self.stat_source)
        else:
            log.info("MST does not need to be revert")
        self._prev_mst = None
        super(MSTchangingMover, self).revert(sm)

class ExhaustiveMover(Mover):
    HELPTEXT = ("Try all fragment combinations for \n"
                "the given coarse-grained element(s).")
    def __init__(self, *elems_of_interest, **kwargs):
        """
        Explore the conformation space exhaustively by iterating
        through the stats for the given elements.

        :param stat_source: A StatContainer instance
        :param sm: The spatial model that will be explored
        :param elems_of_interest: A list of cg-element names.
                        The list gives the order, in which the structure space
                        is explored in a bredth first search.
                        Example: "m0,i1,i2" Try all values for m0, then change
                        i1 once and again try all values of m0, etc
                        Alternatively, a comma-seperated string of element names is supported.

        ..warning ::

            If fallback-stats are used (and there are more fallback-stats than needed),
            the search is not really exhaustive but samples from the fallback stats.
        """
        if "stat_source" not in kwargs:
            raise TypeError("stat_source is a required keyword argument.")
        if "sm" not in kwargs:
            raise TypeError("sm is a required keyword argument")
        stat_source = kwargs["stat_source"]
        sm = kwargs["sm"]
        super(ExhaustiveMover, self).__init__(stat_source)
        if isinstance(elems_of_interest, str):
            elems_of_interest = elems_of_interest.split(",")
        self._element_names = elems_of_interest
        self._move_iterator = self._iter_moves(sm)

    def _iter_moves(self, sm):
        iterables = {}
        for elem in self._element_names:
            iterables[elem]=self.stat_source.iterate_stats_for(sm.bg, elem)
        while True:
            if not iterables:
                return
            for elem in list(iterables.keys()):
                try:
                    stat = next(iterables[elem])
                except StopIteration:
                    del iterables[elem]
                else:
                    yield elem, stat
    def _get_elem_and_stat(self, sm):
        return next(self._move_iterator)



class NElementMover(Mover):
    HELPTEXT = ("In every move, randomly replace \n"
                "n different fragments.")

    def __init__(self, n=2, stat_source=None, **kwargs):
        """
        Change more than one element per iteration step.

        :param n: INT, How many elements should be moved in each step.
        """
        super(NElementMover, self).__init__(stat_source)
        self._n_moves = int(n)
    def move(self, sm):
        if self._n_moves>len(sm.bg.defines):
            raise ValueError("CG RNA has fewer elements than should be change per move.")
        self._prev_stats = {}
        movestring = []
        i=0
        while len(self._prev_stats)<self._n_moves:
            i+=1
            elem, new_stat = self._get_elem_and_stat(sm)
            if elem not in self._prev_stats:
                movestring.append(self._move(sm, elem, new_stat))
            if i>100000:
                raise RuntimeError("Caught in an endless loop (?)")
        sm.new_traverse_and_build(start = "start", include_start=True)
        return "".join(movestring)


class OneOrMoreElementMover(NElementMover):
    HELPTEXT = ("In every move, randomly replace \n"
                "1 to n fragments.")
    def __init__(self, max_n=2, stat_source=None, **kwargs):
        """
        Change more than one element per iteration step.

        :param n: INT, How many elements should be moved in each step.
        """
        super(OneOrMoreElementMover, self).__init__(stat_source)
        self._possible_n = list(range(1, int(max_n)+1))
    @property
    def _n_moves(self):
        return random.choice(self._possible_n)
    @_n_moves.setter
    def _n_moves(self, val):
        # The setter is called in the init function of the superclass.
        # Just ignore it
        pass


class ConnectedElementMover(NElementMover):
    """
    Change more than 1 connected elements per move step
    """
    HELPTEXT = ("In every move, randomly replace \n"
                "n adjacent fragments.")

    def _get_elem(self, sm):
        if not self._prev_stats:
            return super(ConnectedElementMover, self)._get_elem(sm)
        neighbors = { d for changed_elem in self._prev_stats for d in sm.bg.edges[changed_elem]
                        if d not in self._prev_stats
                            and d not in sm.frozen_elements }
        return random.choice(list(neighbors))


    def _get_elem_and_stat(self, sm):
        elem = self._get_elem(sm)
        new_stat = self.stat_source.sample_for(sm.bg, elem)
        return elem, new_stat


class WholeMLMover(Mover):
    """
    Pick a whole multiloop or a random other element and change the stats for it.

    Changing a whole multiloops means changing more than one stat.
    The stats for all ml-segments of this multiloop are sampled independently.
    """
    HELPTEXT = ("In every move, either replace \n"
                "a single not-multiloop fragment,\n"
                "or all multiloop segments \n"
                "for a single multiloop.")

    def move(self, sm):
        movestring = []
        self._prev_stats = {}
        while self._has_incomplete_ml(sm):
            elem, new_stat = self._get_elem_and_stat(sm)
            movestring.append(self._move(sm, elem, new_stat))
        sm.new_traverse_and_build(start = "start", include_start = True)
        return "".join(movestring)
    def _has_incomplete_ml(self, sm):
        if not self._prev_stats:
            return True
        if any(key[0]=="m" for key in self._prev_stats):
            log.info("_has_incomplete_ml: m in prev_stats {}".format(list(self._prev_stats.keys())))
            return self._get_missing_ml_nodes(sm)
        else:
            return False
    def _get_missing_ml_nodes(self, sm):
        sampled = set(k for k in self._prev_stats)
        whole_ml = set(sm.bg.shortest_mlonly_multiloop(list(sampled)[0]))
        missing = whole_ml - sampled - sm.frozen_elements
        log.info("_get_missing_ml_nodes: sampled {},  missing {}".format(sampled, missing))
        return missing

    def _get_elem_and_stat(self, sm):
        if not self._prev_stats:
            return super(WholeMLMover,self)._get_elem_and_stat(sm)
        else:
            elem = random.choice(list(self._get_missing_ml_nodes(sm)))
            new_stat = self.stat_source.sample_for(sm.bg, elem)
            return elem, new_stat
