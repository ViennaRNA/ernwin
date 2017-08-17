#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input, #pip install future
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)
__metaclass__=type
import itertools as it
import random
import logging
import inspect
import math

try:
    from types import SimpleNamespace
except ImportError: #<python3.3
    class SimpleNamespace:
        pass

import numpy as np

import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.vector as ftuv
import forgi.threedee.model.stats as ftms

from ..aux.utils import get_all_subclasses

log=logging.getLogger(__name__)

class Mover:
    HELPTEXT = ("{:25} Randomly replace a single fragment\n"
                "{:25} of the RNA.".format("Mover", ""))
    def __init__(self, stat_source):
        """
        A Mover class that randomly pickes an element and changes the stats for it.
        """
        self.stat_source = stat_source
        #: A list of tuples  (elemenmt_name, stat)
        self._prev_stats = None

    def _get_elem(self, sm):
        possible_elements = sm.bg.get_mst() - sm.frozen_elements
        return random.choice(list(possible_elements))

    def _get_elem_and_stat(self, sm):
        elem = self._get_elem(sm)
        new_stat = self.stat_source.sample_for(sm.bg, elem)
        return elem, new_stat

    def move(self, sm):
        """
        Propose a single Monte Carlo move.

        :param sm: The spatial model
        """
        elem, new_stat = self._get_elem_and_stat(sm)
        self._prev_stats = {}
        movestring = self._move(sm, elem, new_stat)
        sm.new_traverse_and_build(start = elem, include_start = True)
        return movestring

    def _move(self, sm, elem, new_stat):
        try:
            prev_stat = sm.elem_defs[elem]
        except KeyError:
            prev_stat = SimpleNamespace()
            prev_stat.pdb_name = "UNSET"
        else:
            self._prev_stats[elem] = prev_stat
        sm.elem_defs[elem]=new_stat
        return "{}:{}->{};".format(elem, prev_stat.pdb_name, new_stat.pdb_name)

    def revert(self, sm):
        """
        Revert the last Monte Carlo move performed by this mover.
        """
        log.debug("Reverting last step")
        if self._prev_stats is None:
            raise RuntimeError("No (more) step(s) to revert.")
        for elem, stat in self._prev_stats.items():
            sm.elem_defs[elem] = stat
        self._prev_stats = {}
        sm.new_traverse_and_build(start='start', include_start = True)

class MSTchangingMover(Mover):
    def __init__(self, stat_source):
        super(self, MSTchangingMover).__init__(stat_source)
        self.prev_mst = None
    def _get_elem_and_stat(self, sm):
        # Get an element. If it is a ml-segment, break it.
        elem = self._get_elem()
        if elem[0]=="m":
            elem = sm.set_multiloop_break_segment(elem)
        new_stat = self.stat_source.sample_for(sm.bg, elem)
        return elem, new_stat



class ExhaustiveMover(Mover):
    HELPTEXT = ("{:25} Try all fragment combinations for \n"
                "{:25} the given coarse-grained element(s).".format("ExhaustiveMover(elem[,elem,..])", ""))
    def __init__(self, stat_source, elems_of_interest, sm):
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
        super(ExhaustiveMover, self).__init__(stat_source)
        if isinstance(elems_of_interest, str):
            elems_of_interest = elems_of_interest.split(",")
        print(elems_of_interest)
        self._element_names = elems_of_interest
        self._move_iterator = self._iter_moves(sm)

    def _iter_moves(self, sm):
        iterables = []
        for elem in self._element_names:
            iterables.append(zip(it.repeat(elem),
                                 self.stat_source.iterate_stats_for(sm.bg, elem)))
        old_stat  = {}
        for prod in it.product(*iterables):
            for elem, stat in prod:
                if elem in old_stat and old_stat[elem] == stat:
                    continue
                else:
                    old_stat[elem]=stat
                    yield elem, stat

    def _get_elem_and_stat(self, sm):
        return next(self._move_iterator)

class NElementMover(Mover):
    HELPTEXT = ("{:25} In every move, randomly replace \n"
                "{:25} n different fragments.".format("NElementMover(N)", ""))

    def __init__(self, stat_source, n=2):
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
    HELPTEXT = ("{:25} In every move, randomly replace \n"
                "{:25} 1 to n fragments.".format("OneOrMoreElementMover[N]", ""))
    def __init__(self, stat_source, max_n=2):
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
    HELPTEXT = ("{:25} In every move, randomly replace \n"
                "{:25} n adjacent fragments.".format("ConnectedElementMover(N)", ""))

    def _get_elem(self, sm):
        if not self._prev_stats:
            return super(ConnectedElementMover, self)._get_elem(sm)
        neighbors = { d for changed_elem in self._prev_stats for d in sm.bg.edges[changed_elem]
                        if d not in self._prev_stats
                            and d in sm.bg.mst
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
    HELPTEXT = ("{:25} In every move, either replace \n"
                "{:25} a single not-multiloop fragment,\n"
                "{:25} or all multiloop segments \n"
                "{:25} for a single multiloop.".format("WholeMLMover", "", "", ""))

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
            log.info("_has_incomplete_ml: m in prev_stats {}".format(self._prev_stats.keys()))
            return self._get_missing_ml_nodes(sm)
        else:
            return False
    def _get_missing_ml_nodes(self, sm):
        sampled = set(k for k in self._prev_stats)
        whole_ml = set(sm.bg.shortest_mlonly_multiloop(list(sampled)[0]))
        whole_ml = whole_ml & sm.bg.mst
        missing = whole_ml - sampled - sm.frozen_elements
        log.info("_get_missing_ml_nodes: sampled {}, mst {}, missing {}".format(sampled, sm.bg.mst, missing))
        return missing

    def _get_elem_and_stat(self, sm):
        if not self._prev_stats:
            return super(WholeMLMover,self)._get_elem_and_stat(sm)
        else:
            elem = random.choice(list(self._get_missing_ml_nodes(sm)))
            new_stat = self.stat_source.sample_for(sm.bg, elem)
            return elem, new_stat

class WholeMLStatSearch(Mover):
    HELPTEXT = ("{:25} Not yet implemented\n".format("WholeMLStatSearch"))

    def __init__(self, stat_source, max_diff):
        super(WholeMLStatSearch, self).__init__(stat_source)
        self.max_diff=float(max_diff)
    def _get_elems(self, sm):
        loops = sm.bg.find_mlonly_multiloops()
        regular_multiloops = [ m for m in loops
                               if "regular_multiloop" in sm.bg.describe_multiloop(m) ]
        if len(regular_multiloops)==0:
            raise ValueError("{} needs at least 1 regular multiloop. "
                             "(Pseudoknots and external loops are not yet supported.)")
        return random.choice(regular_multiloops)
    def move(self, sm):
        self.counter = 0
        self.num_pruned = 0
        elems = self._get_elems(sm)
        assert sm.bg.get_next_ml_segment(elems[0])==elems[1]
        assert sm.bg.get_next_ml_segment(elems[1])==elems[2]
        assert sm.bg.get_next_ml_segment(elems[-1])==elems[0]

        stats = self._find_stats_for(elems, sm)
        if not stats:
            log.warning("Could not move multiloop %s. "
                        "No suitable combination of stats found.", elems)
            return "no_change"
        else:
            movestring = []
            self._prev_stats = {}
            for ml, stat in stats.items():
                movestring.append(self._move(sm, ml, stat))
            sm.new_traverse_and_build(start = "start", include_start = True)
            return "".join(movestring)
    def _find_stats_for(self, elems, sm):
        """
        :param elems: ML segments. Need to be ordered!
        """
        choices = { elem : list(self.stat_source.iterate_stats_for(sm.bg, elem)) for elem in elems }
        for stat_choices in choices.values():
            random.shuffle(stat_choices) # For performance reasons lists are shuffled in place
        maxlen = max(len(c) for c in choices.values())
        log.debug("maxlen %s", maxlen)
        for i in range(maxlen):
            for first_elem in choices.keys():
                sampled_stats = self._find_stats_ith_iteration(i, first_elem, elems, choices, sm.bg)
                if sampled_stats is not None:
                    log.info("Succsessfuly combination found after %d tries in %d iterations (%d of which were pruned)", self.counter, i, self.num_pruned)
                    return sampled_stats
            if i%10==9:
                log.info("Nothing found after %d iterations. %d combinations tried, %d of which were disregarded early on. Still searching.", i+1, self.counter, self.num_pruned)
        return None
    def _find_stats_ith_iteration(self, i, first_elem, elems, choices, cg):
        first_elem_index = elems.index(first_elem)
        # Of the selected ml-segment, only look at the ith stat and compare it
        # to all possible combinations of the first i stats for other elements.
        try:
            first_elem_stat = choices[first_elem][i]
        except IndexError:
            return None
        log.debug("i = %d. elems = %s, first_elem %s", i, elems, first_elem)
        for choice_indices in it.product(range(i+1), repeat=len(elems)-1):
            log.debug("Choice_indices %s", choice_indices)
            sampled_stats = {}
            # For the ml-segments AFTER the first_elem, only try stats until the (i-1)th
            log.debug("Range (%d, %d): %s", first_elem_index, len(elems)-1, list(range(first_elem_index, len(elems)-1)))
            if any(choice_indices[elem_nr]==i for elem_nr in range(first_elem_index, len(elems)-1)):
                log.debug("continue")
                continue
            skip = False
            for elem_nr in range(len(choice_indices)):
                if elem_nr>=first_elem_index:
                    log.debug("elem_nr == %d --> %d", elem_nr, elem_nr+1)
                    ml = elems[elem_nr+1]
                else:
                    log.debug("elem_nr %d", elem_nr)
                    ml = elems[elem_nr]
                assert ml != first_elem
                try:
                    new_stat = choices[ml][choice_indices[elem_nr]]
                except IndexError:
                    skip=True
                    break
                sampled_stats[ml] = new_stat
            if not skip:
                sampled_stats[first_elem] = first_elem_stat
                if self._check_overall_stat(sampled_stats, elems, cg):
                    return sampled_stats
        return None
    def _check_overall_stat(self, sampled_stats, elems, cg):
        self.counter+=1
        log.debug("Checking sampled_stats %s for elems %s", sampled_stats, elems)
        partial_sum_stat = None
        if len(elems)==3:
            #speed-up by pruning according to triangular inequality
            dists = [stat.r1 for stat in sampled_stats.values()]
            sum_of_lengths = sum(dists)
            max_length = max(dists)
            if sum_of_lengths - max_length> max_length+self.max_diff:
                log.debug("Pruning because of triangular inequality: sum %s, max %s", sum_of_lengths, max_length )
                self.num_pruned+=1
                return False
        for elem in elems:
            at = cg.get_angle_type(elem, allow_broken = True)
            if at in [-3, -2, 4]:
                next_stat = - sampled_stats[elem]
            else:
                next_stat = sampled_stats[elem]
            if partial_sum_stat is None:
                partial_sum_stat = next_stat
            else:
                partial_sum_stat = ftug.sum_of_stats(partial_sum_stat, next_stat)

        return partial_sum_stat.is_similar_to(ftms.IDENTITY_STAT, self.max_diff, math.radians(self.max_diff))

class LocalMLMover(Mover):
    """
    Change two connected ML elements in a way that does not change the parts
    of the RNA attached to the two ends of the connected elements.
    """
    HELPTEXT = ("{:25} Not yet implemented\n".format("LocalMLMover"))

    def __init__(self, stat_source, max_diff):
        super(LocalMLMover, self).__init__(stat_source)
        self.max_diff=float(max_diff)

    def _get_elems(self, sm):
        all_mls = list( sm.bg.mloop_iterator() )
        if len(all_mls)<2:
            raise ValueError("LocalMLMover needs at least 2 multiloop"
                             " segments in the sm")
        ml1 = random.choice(all_mls)
        ml2 = sm.bg.get_next_ml_segment(ml1)
        log.info("Moving elements %s and %s. MST is %s.", ml1, ml2, sm.bg.mst)
        log.info("Elem_defs for: %s", list(sm.elem_defs.keys()))
        return ml1, ml2

    def _find_stats_ith_iteration(self, i, choices1, choices2, virtual_stat, forward= True):
        """
        Compare the ith stat of the first choices,
        to the first until ith stat of the second choices
        """
        if forward:
            choicesA, choicesB= choices1, choices2
        else:
            choicesB, choicesA= choices1, choices2

        if i>=len(choicesA):
            return None
        statA = choicesA[i]
        for j in range(min(i+forward, len(choicesB))):
            statB = choicesB[j]
            # Shortcut based on lengths
            if virtual_stat.r1>statA.r1+statB.r1+self.max_diff:
                return None
            if forward:
                if virtual_stat.is_similar_to(ftug.sum_of_stat_in_standard_direction(statA, statB), self.max_diff):
                    return statA, statB
            else:
                if virtual_stat.is_similar_to(ftug.sum_of_stat_in_standard_direction(statB, statA), self.max_diff):
                    return statB, statA
        return None

    def _find_stats_for(self, ml1, ml2, sm):
        virtual_stat = ftug.get_virtual_stat(sm.bg, ml1, ml2)

        stem1a, stem1b = sm.bg.connections(ml1)
        stem2a, stem2b = sm.bg.connections(ml2)
        swapped = False
        if stem1b == stem2b:
            swapped = True
            ml1, ml2 = ml2, ml1

        log.debug("Finding stats for %s and %s", ml1, ml2)
        choices1 = list(self.stat_source.iterate_stats_for(sm.bg, ml1))
        choices2 = list(self.stat_source.iterate_stats_for(sm.bg, ml2))
        random.shuffle(choices1)
        random.shuffle(choices2)
        maxlen = max(len(choices1), len(choices2))
        try:
            for i in range(maxlen):
                for forward in [True, False]:
                    result = self._find_stats_ith_iteration(i, choices1, choices2, virtual_stat, forward)
                    if result is not None:
                        log.info("Stat found after %sth iteration. Forward==%s", i, forward)
                        if swapped:
                            log.debug("Swapping...")
                            return result[1], result[0]
                        else:
                            return result
                if i>0 and i%50==0:
                    log.info("Search of stats: Nothing found in %d th iteration", i)

        except Exception as e:
            log.error("Interrupted in %sth iteration (out of %s)", i, maxlen)
            raise
        return None, None

    def move(self, sm):
        ml1, ml2 = self._get_elems(sm)

        if sm.bg.define_a(ml2)[0]<sm.bg.define_a(ml1)[0]:
            ml1, ml2 = ml2, ml1

        stat1, stat2 = self._find_stats_for(ml1, ml2, sm)
        log.debug("%s %s; %s %s", ml1, stat1, ml2, stat2)
        if stat1 is stat2 is None:
            log.warning("Could not move elements %s and %s. "
                        "No isosteric combination of stats found.", ml1, ml2)
            return "no_change"
        else:
            movestring = []
            self._prev_stats = {}
            log.info("Moving %s, elem_defs=%s", ml1, list(sm.elem_defs.keys()))
            movestring.append(self._move(sm, ml1, stat1))
            log.info("Moving %s, elem_defs=%s", ml2, list(sm.elem_defs.keys()))
            movestring.append(self._move(sm, ml2, stat2))
            sm.new_traverse_and_build(start = "start", include_start = True)
            return "".join(movestring)


####################################################################################################
### Command line parsing
####################################################################################################

def mover_from_string(stri, stat_source, sm=None):
    """
    Return a single Mover instance from a string describing the mover.

    The string needs to contain the class name of the mover and optionally
    may contain one argument in square brackets: MOVERNAME[ARGUMENT]

    :param stat_source: The stat_container that shall be used for all moves.
    :param sm: The Spatiel model. Needed, for ExhaustiveMover

    """
    if '[' in stri:
        movername, option = stri.split('[')
        if option[-1]!="]":
            raise ValueError("Missing closing bracket at end of '[{}'".format(option))
        args = [option[:-1]]
    else:
        args = []
        movername = stri
    kwargs = {}
    mover_classes = { cls.__name__: cls for cls in get_all_subclasses(Mover, include_base = True) }
    try:
        cls = mover_classes[movername]
        required = inspect.getargspec(cls.__init__).args
        if "sm" in required:
            kwargs["sm"] = sm
        mover = cls(stat_source, *args, **kwargs)
    except TypeError as e: #Give more helpful error messages
        if "arguments" not in e.message or "__init__" not in e.message:
                raise
        # Something like TypeError: __init__() takes exactly 4 arguments (3 given)
        signature = inspect.getargspec(cls.__init__)
        called = "[self], stat_source,"
        called+=",".join(map(str,args))
        if kwargs:
            called+=","
            for k,v in kwargs.items():
                called+="{}={},".format(k,v)
        raise TypeError("__init__ of {} with signature ({}) cannot be called with ({})".format(cls.__name__, signature, called))

    return mover

def get_argparse_help():
    """
    Return a pre-formatted string that can be passed to "help" in argparse.add_argument,
    if the resulting argument is parsed with `mover_from_string`

    .. warning::

        Please specify the same default in `add_argument` as described in this help-text
    """
    help= ("Which types of Moves to use during sampling.\n"
           "Default: Mover\n"
           "One of the following:\n")
    for mover_class in get_all_subclasses(Mover, include_base = True):
        help+=mover_class.HELPTEXT+"\n"
    return help
