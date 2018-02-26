#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input, #pip install future
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
__metaclass__=type
import itertools as it
import random
import logging
import inspect
import math
import copy
import textwrap
from commandline_parsable import parsable_base

try:
    from types import SimpleNamespace
except ImportError: #<python3.3
    class SimpleNamespace:
        pass

import numpy as np

import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.vector as ftuv
import forgi.threedee.model.stats as ftms

from ..utils import get_all_subclasses
from . import create
from ._commandline_helper import replica_substring


log=logging.getLogger(__name__)

class UnsuitableMover(ValueError):
    """
    Raised by a Mover if it is called by an RNA which does not have
    the SECONDARY structure required for this mover to work.

    For example, a Mover that works on multiloops will raise this,
    if the structure has no multiloop.

    .. note:: This is NOT raised, if the mover fails to move
              the structure due tertiary structure clashes, because
              these clashes might be removed by first moving other
              parts of the structure.
    """

@parsable_base(True, required_kwargs=["stat_source"], helptext_sep="\n",
               help_attr="HELPTEXT", help_intro_list_sep="\n")
class Mover:
    HELPTEXT = ("Randomly replace a single fragment\n"
                "of the RNA.")
    def __init__(self, stat_source, **kwargs):
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
        log.info("%s move called", type(self).__name__)
        elem, new_stat = self._get_elem_and_stat(sm)
        self._prev_stats = {}
        movestring = self._move(sm, elem, new_stat)
        sm.new_traverse_and_build(start = elem, include_start = True)
        return movestring

    def _store_prev_stat(self, sm, elem):
        """
        Store the stat from elem_defs to self._prev_stats and return its pdb_name.
        Returns "UNSET" if the stat is not present.
        """
        try:
            prev_stat = sm.elem_defs[elem]
        except KeyError:
            return "UNSET"
        else:
            self._prev_stats[elem] = prev_stat
        return prev_stat.pdb_name
    def _move(self, sm, elem, new_stat):
        prev_name = self._store_prev_stat(sm, elem)
        sm.elem_defs[elem]=new_stat
        return "{}:{}->{};".format(elem, prev_name, new_stat.pdb_name)

    def revert(self, sm):
        """
        Revert the last Monte Carlo move performed by this mover.
        """
        log.debug("%s Reverting last step", type(self).__name__)
        if self._prev_stats is None:
            raise RuntimeError("No (more) step(s) to revert.")
        for elem, stat in self._prev_stats.items():
            assert stat is not None
            sm.elem_defs[elem] = stat
        self._prev_stats = {}
        sm.new_traverse_and_build(start='start', include_start = True)


class MoverNoRegularML(Mover):
    def _get_elem(self, sm):
        while True:
            elem = super(MoverNoRegularML, self)._get_elem(sm)
            if elem[0] == "m":
                loop = sm.bg.shortest_mlonly_multiloop(elem)
                if "regular_multiloop" not in sm.bg.describe_multiloop(loop):
                    return elem
            else:
                return elem



class EnergeticJunctionMover(Mover):
    HELPTEXT= textwrap.fill(textwrap.dedent("""\
                Try sampling for n consecutive fragments of a junction,
                until the junction constraint energy (if any)
                for this junction is fulfilled. Use n=-1 for
                the whole junction"""), 35)
    def __init__(self, n, stat_source, **kwargs):
        if n!=-1 and n<1:
            raise ValueError("EnergeticJunctionMover requires n==-1 or n>=1, not n={}".format(n))
        super(EnergeticJunctionMover, self).__init__(stat_source)
        self.n_elements = n
        self.max_tries = 20000
        self.original_max_tries = self.max_tries
        #: The first time move is called, we enumerate all choices,
        #: so we can rule out the possibility of 0 available choices.
        #: A dictionary sm-identifier : choices.
        self.choices = {}
    @staticmethod
    def _sm_fingerprint( sm):
        fingerprint = sm.bg.name+"\n"+sm.bg.to_dotbracket_string()+"\n"
        fingerprint+=",".join(sorted(sm.frozen_elements))
        return fingerprint
    def _enumerate_choices(self, sm):
        choices = []
        loops = sm.bg.find_mlonly_multiloops()
        regular_multiloops = [ m for m in loops
                               if "regular_multiloop" in sm.bg.describe_multiloop(m) ]
        if len(regular_multiloops)==0:
            raise UnsuitableMover("{} needs at least 1 regular multiloop. "
                                  "(Pseudoknots and external loops are "
                                  "not yet supported.)")
        for loop in regular_multiloops:
            if self.n_elements==-1:
                if all(m not in sm.frozen_elements for m in loop):
                    choices.append(self._sort_loop(sm, list(loop)))
            elif self.n_elements>0:
                for elem in loop:
                    elements = [elem]
                    for i in range(1, self.n_elements):
                        next_elem = sm.bg.get_next_ml_segment(elements[-1])
                        if next_elem in elements:
                            break
                        else:
                            elements.append(next_elem)
                    if elements not in choices and all(e not in sm.frozen_elements
                                                       for e in elements):
                        choices.append(elements)
            else:
                assert False
        return choices

    def _get_elements(self, sm):
        fingerprint = self._sm_fingerprint(sm)
        if fingerprint not in self.choices:
            c = self._enumerate_choices(sm)
            if len(c)==0:
                raise UnsuitableMover("Cannot move any regulkr multiloop, "
                                      "because all contain frozen elements.")
            self.choices[fingerprint] = c
        return random.choice(self.choices[fingerprint])

    def _sort_loop(self, sm, loop):
        """
        Sort loop according to buildorder
        """
        build_order = sm.bg.traverse_graph()
        def s_key(ml):
            for i, stem_loop_stem in enumerate(build_order):
                if ml==stem_loop_stem[1]:
                    return i
            return float("inf")
        loop.sort(key=s_key)
        log.debug("Loop sorted: %s", loop)
        return loop

    def move(self, sm):
        elems = self._get_elements(sm)

        self._prev_stats = {}
        for elem in elems:
            self._store_prev_stat(sm, elem)

        use_asserts = ftuv.USE_ASSERTS
        ftuv.USE_ASSERTS=False
        stats = self._find_stats_for(elems, sm)
        ftuv.USE_ASSERTS=use_asserts

        if not stats:
            log.warning("Could not move multiloop %s. "
                        "No suitable combination of stats found.", elems)
            self.revert(sm)
            return "no_change"
        else:
            movestring = []

            for ml, stat in stats.items():
                try:
                    p_name = self._prev_stats[ml].pdb_name
                except KeyError:
                    p_name = "UNSET"
                movestring.append("{}:{}->{};".format(ml, p_name, stat.pdb_name))
            sm.new_traverse_and_build(start = "start", include_start = True)
            return "".join(movestring)

    def _check_junction(self, sm, sampled_stats, elems):
        for elem, new_stat in sampled_stats.items():
            sm.elem_defs[elem]=new_stat
        built_nodes = []
        for elem in elems[:-1]: # The last elem is broken!
            nodes = sm.new_traverse_and_build(start=elem, max_steps = 1, finish_building=False)
            built_nodes+=nodes
        try:
            broken_stat = sampled_stats[elems[-1]]
        except KeyError:
            try:
                broken_stat = sm.elem_defs[elems[-1]]
            except KeyError: # We use the JDIST energy
                broken_stat = None
        energy = sm.junction_constraint_energy[elems[-1]].eval_energy(sm.bg,
                                                            nodes=built_nodes+[elems[-1]],
                                                            sampled_stats={elems[-1]:broken_stat})
        if energy==0:
            return True
        return False

    def _find_stats_for(self, elems, sm):
        """
        :param elems: ML segments. Need to be ordered!
        """
        counter = 0
        loop = self._sort_loop(sm, list(sm.bg.shortest_mlonly_multiloop(elems[0])))
        # We do not have to build the whole loop,
        # but only the part after the first changed element.
        i = min(loop.index(elem) for elem in elems)
        log.info("For elems %s, loop index is %s, loop is %s", elems, i, loop)
        loop = loop[i:]
        log.info("Loop now %s", loop)
        for sampled in create.stat_combinations(sm.bg, elems, self.stat_source):
            counter+=1
            if self._check_junction(sm, sampled, loop):
                log.info("Succsessfuly combination found after %d tried", counter)
                if self.max_tries is not None:
                    self.max_tries = int(max(self.original_max_tries, self.max_tries*3/4, counter+(self.original_max_tries/2)))
                    log.info("Setting self.max_tries to %d", self.max_tries)
                return sampled
            if counter%10000==0 and counter>0:
                log.info("Nothing found after %d tries for %s."
                         "Still searching.", counter, elems)
            if counter>self.max_tries:
                # We give up without having exhausted the search space, but we will
                # try twice as hard next time.
                # Since we might have different ml-segments next time,
                # it may be faster anyway.
                # This way we avoid spending unproportionally long time in one
                # step, if all other steps would be quick.
                log.info("Giving up after %d tries this time, setting max_tries to %s.", counter, self.max_tries*2)
                self.max_tries*=2
                return None
        return None

class MixedMover():
    def __init__(self, movers=[]):
        self.movers = movers
        self.last_mover = None
    def move(self, sm):
        try:
            self.last_mover = random.choice(self.movers)
        except IndexError:
            raise ValueError("None of the provided Movers is suitable for this RNA.")
        try:
            return self.last_mover.move(sm)
        except UnsuitableMover:
            self.movers.remove(self.last_mover)
            return self.move(sm)

    def revert(self, sm):
        self.last_mover.revert(sm)


####################################################################################################
### Command line parsing
####################################################################################################

def _mover_from_string(stri, stat_source, sm=None):
    # Use the from_string classmethod provided by @parsable_base
    log.debug("Generating Mover from string %r, with stat_source=%r, sm=%r", stri, stat_source, sm)
    movers = Mover.from_string(stri, stat_source=stat_source, sm=sm)
    if len(movers)==0:
        raise ValueError("At least one Mover has to be specified.")
    elif len(movers)==1:
        return movers[0]
    else:
        return MixedMover(movers)

def update_parser(parser):
    mover_help_intro = ("Which types of Moves to use during sampling.\n"
                        "If you specify more than one mover separated by a\n"
                        "comma, the mover will be picked at random at each step.\n"
                        "Default: Mover\n"
                        "One or more of the following:")
    Mover.add_to_parser(parser, '--move-set', default="MoverNoRegularML,EnergeticJunctionMover[2]", help_intro=mover_help_intro)

def from_args(args, stat_source, sm=None, replica_nr = None):
    argument_string = replica_substring(args.move_set, replica_nr)
    return _mover_from_string(argument_string, stat_source, sm)
