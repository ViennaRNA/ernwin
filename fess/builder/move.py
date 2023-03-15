#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input, #pip install future
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from six.moves import map
from six.moves import range
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
from . import relaxation_builder as fbrel

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
        possible_elements = set(sm.bg.defines.keys()) - sm.frozen_elements
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
        if elem in sm.bg.get_mst():
            start=elem
        else:
            start="end"
        log.debug("Mover building (start=%s")
        sm.new_traverse_and_build(start = start, include_start = True)
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
        log.debug("%s MOVE: Assigning %s to %s", type(self).__name__, new_stat.pdb_name, elem)
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
            log.debug("%s REVERT Assigning %s to %s", type(self).__name__, stat.pdb_name, elem)
            sm.elem_defs[elem] = stat
        self._prev_stats = {}
        sm.new_traverse_and_build(start='start', include_start = True)

class MoveAndRelaxer(Mover):
    def _store_prev_stat(self, sm, elem):
        """
        Store the stat from elem_defs to self._prev_stats and return its pdb_name.
        Returns "UNSET" if the stat is not present.

        Due to relaxation potentially affecting all stats, we just store all in this mover.
        """
        for e, stat in sm.elem_defs.items():
            if e != elem:
                self._prev_stats[e] = stat
        try:
            prev_stat = sm.elem_defs[elem]
        except KeyError:
            return "UNSET"
        else:
            self._prev_stats[elem] = prev_stat
        return prev_stat.pdb_name

    def move(self, sm):
        """
        Propose a single Monte Carlo move.
        And relax the rest of the structure to remove clashes.

        :param sm: The spatial model
        """
        log.info("%s move called", type(self).__name__)
        elem, new_stat = self._get_elem_and_stat(sm)
        self._prev_stats = {}
        movestring = self._move(sm, elem, new_stat)
        if elem in sm.bg.get_mst():
            start=elem
        else:
            start="end"
        sm.new_traverse_and_build(start = start, include_start = True)
        ok, relaxstring = fbrel.relax_sm(sm, self.stat_source, [elem])
        return movestring+relaxstring

class MoveAndRelaxML(MoveAndRelaxer):
    def _get_elem(self, sm):
        possible_elements = set(sm.bg.mloop_iterator()) - sm.frozen_elements
        if not possible_elements:
            raise UnsuitableMover("No ML in structure")
        return random.choice(list(possible_elements))

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
        n=int(n)
        if n!=-1 and n<1:
            raise ValueError("EnergeticJunctionMover requires n==-1 or n>=1, not n={}".format(n))
        super(EnergeticJunctionMover, self).__init__(stat_source)
        self.n_elements = n
        if n==-1:
            self.max_tries = None
        else:
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
            log.debug("Enumerating choices for loop %s with lengths %s",
                        loop, list(map(sm.bg.element_length, loop)))
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
                raise UnsuitableMover("Cannot move any regular multiloop, "
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
        return self.move_elems(sm, elems)

    def move_elems(self, sm, elems):
        """
        Called directly by DimerizationBuilder.
        """
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

    def _check_junction_pk():
        """
        No clever optimization for Pseudoknots. Just build the whole structure.
        """
        # TODO Implement

    def _check_junction(self, sm, sampled_stats, elems, whole_loop):
        for elem, new_stat in sampled_stats.items():
            sm.elem_defs[elem]=new_stat
        # Energy may implement a pre-check to rule-out impossible conformations
        # This heuristic saves 10-15% of evaluations for some text structures
        # and saves between 0 and 10% of the runtime ...
        try:
            if sm.junction_constraint_energy[elems[-1]].precheck(sm.bg, whole_loop, sm.elem_defs)>0:
                return False
        except AttributeError:
            pass
        # Now we have to build the junction
        built_nodes = []
        for elem in elems[:-1]: # The last elem is broken!
            nodes = sm.new_traverse_and_build(start=elem, max_steps = 1)#, finish_building=False)
            built_nodes+=nodes
        try:
            broken_stat = sampled_stats[elems[-1]]
        except KeyError:
            try:
                broken_stat = sm.elem_defs[elems[-1]]
            except KeyError: # We use the JDIST energy
                broken_stat = None
        log.debug("sm has junction-energy for %s. Elem is %s",
                  list(sm.junction_constraint_energy.keys()), elems[-1])
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
        whole_loop = self._sort_loop(sm, list(sm.bg.shortest_mlonly_multiloop(elems[0])))
        # We do not have to build the whole loop,
        # but only the part after the first changed element.
        i = min(whole_loop.index(elem) for elem in elems)
        log.debug("For elems %s, loop index is %s, loop is %s", elems, i, whole_loop)
        loop = whole_loop[i:]
        log.info("Loop now %s", loop)
        for sampled in create.stat_combinations(sm.bg, elems, self.stat_source):
            counter+=1
            if self._check_junction(sm, sampled, loop, whole_loop):
                log.info("Succsessfuly combination found after %d tried", counter)
                if self.max_tries is not None:
                    self.max_tries = int(max(self.original_max_tries, self.max_tries*3/4, counter+(self.original_max_tries/2)))
                    log.info("Setting self.max_tries to %d", self.max_tries)
                return sampled
            if counter%10000==0 and counter>0:
                log.info("Nothing found after %d tries for %s."
                         "Still searching.", counter, elems)
            if self.max_tries is not None and counter>self.max_tries:
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


class RotationMover(Mover):
    def __init__(self, stat_source, **kwargs):
        self.last_axis = None
        self.last_angle = None

    def move(self, sm):
        self.last_axis = random.choice(["x", "y", "z"])
        self.last_angle = random.random()*2*np.pi # 0-360 degrees
        rot_mat = ftuv.rotation_matrix(self.last_axis, self.last_angle)
        sm.bg.rotate(rot_mat)
        return "{}deg{}".format(math.degrees(self.last_angle), self.last_axis)

    def revert(self, sm):
        assert self.last_axis is not None
        sm.bg.rotate(self.last_axis, -self.last_angle)
        self.last_axis = None


class MixedMover():
    def __init__(self, movers=[]):
        self.movers = movers
        self.moveAndRelML = None
        for mover in self.movers:
            if isinstance(mover, MoveAndRelaxML):
                self.moveAndRelML = mover
                self.movers.remove(mover)
        self.last_mover = None
        self.i=0
    def move(self, sm):
        self.i+=1
        try:
            if self.i%25==1 and self.moveAndRelML is not None:
                self.last_mover = self.moveAndRelML
            else:
                self.last_mover = random.choice(self.movers)
        except IndexError:
            raise ValueError("None of the provided Movers is suitable for this RNA.")
        try:
            return self.last_mover.move(sm)
        except UnsuitableMover:
            try:
                self.movers.remove(self.last_mover)
            except ValueError:
                self.moveAndRelML = None
            self.i-=1
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
                        "Default: 'MoverNoRegularML,EnergeticJunctionMover[2]'\n"
                        "One or more of the following:")
    Mover.add_to_parser(parser, '--move-set', default="MoverNoRegularML,EnergeticJunctionMover[2]", help_intro=mover_help_intro)

def from_args(args, stat_source, sm=None, replica_nr = None):
    argument_string = replica_substring(args.move_set, replica_nr)
    return _mover_from_string(argument_string, stat_source, sm)
