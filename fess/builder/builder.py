#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)
"""builder.py: This file contains classes, which take a spatial model without any 3D information
and add 3D information to it."""

__author__ = "Bernhard Thiel"
__copyright__ = "Copyright 2016"
__license__ = "GNU Affero GPL v 3.0"
__maintainer__ = "Bernhard Thiel"
__email__ = "thiel@tbi.univie.ac.at"

import itertools
import random
import copy
import os
import sys
#profile decorator from line_profiler (kernprof) or memory_profiler
try:
    profile
except:
    profile = lambda x: x

import forgi.threedee.utilities.graph_pdb as ftug
import logging
log = logging.getLogger(__name__)

def load_sampled_elements(sm):
    """
    Try to load the sampled elements from the cg into the sm.
    :returns: True upon success, False upon failure (e.g. if cg was derived from a fasta file)
    """
    try:
        sm.load_sampled_elems()
    except Exception as e:
        log.warning("Cannot use structure from file. Need to resample.")
        log.debug("Reason for need to resampe: {}".format(e), exc_info=True)
        return False
    if not sm.elem_defs:
        return False
    return True

def _determined_broken_ml_segments(built_nodes, bg):
    """
    Return a list of (broken) ml segments that are fully determined by the nodes
    in built_nodes.
    """
    ml_nodes=set(x for x in bg.defines.keys() if x[0]=="m")
    broken_multiloops = ml_nodes-set(itertools.chain(*[bo for bo in bg.traverse_graph()]))
    log.debug("MST = %s, build_order= %s", bg.mst, bg.build_order)
    log.debug("Broken determined multiloops are %s. Now finding out if they are determined...", broken_multiloops)
    broken_determined_nodes=set()
    for n in broken_multiloops:
        loop=set(bg.find_bulge_loop(n, 200))
        if loop: #loop is empty for ml at 3'/ 5' end.
            if loop <= ( set(built_nodes) | broken_multiloops ):
                broken_determined_nodes.add(n)
    return broken_determined_nodes


class Builder(object):
    """
    Build a structure with arbitrary stats from the stat-source,
    that fulfill all constraint energies but are
    not representative samples from the conformation space.
    """
    def __init__(self, stat_source):
        self.stat_source = stat_source
        self.clash_only_tries = 100

    def build_n(self, sm, n):
        """
        Return a list of initialized copies of the spatial model

        :param sm: The spatial model
        :param n: The number of builds you like to get.
        """
        models = []
        for i in range(n):
            self.build(sm)
            models.append(copy.deepcopy(sm))
        return models

    def accept_or_build(self, sm):
        log.info("building without constraint energy...")
        if not sm.elem_defs:
            return self.build(sm)
        else:
            sm.new_traverse_and_build()
            if sm.constraint_energy is not None or sm.junction_constraint_energy:
                self._build_with_energies(sm)
            log.info("Done to build")

    def build(self, sm):
        """
        Initialize the spatial model with random stats

        :param sm: The SpatialModel
        """
        sm.sample_stats(self.stat_source)
        self.accept_or_build(sm)

    def _build_with_energies(self, sm):
        log.info("building with constraint energies")

        newbuilt_nodes = sm.new_traverse_and_build(start = 'start', max_steps = 1)
        built_nodes = []
        iterations = 0
        try:
            while newbuilt_nodes:
                iterations +=1
                built_nodes += newbuilt_nodes
                log.debug("built nodes are {}".format(built_nodes))
                self._check_sampled_ml(sm, newbuilt_nodes[-1])
                bad_segments = self._get_bad_ml_segments(sm, built_nodes)
                log.debug("Bad segments {}".format(bad_segments))
                if not bad_segments:
                    log.debug("Evaluate clash energy:")
                    #The junction-energy is ok. Now we look at the clash energy
                    bad_segments = self._get_bad_clash_segments(sm, built_nodes)
                    if not bad_segments or self._rebuild_clash_only(sm, built_nodes, [x for x in bad_segments if x[0] == "i"]):
                        log.debug("clashfree.")
                        #All clashes were removed
                        assert self._get_bad_clash_segments(sm, built_nodes) == []
                        bad_segments = []
                    else:
                        #The structure has changed, so we need to get the bad segments again.
                        bad_segments = self._get_bad_clash_segments(sm, built_nodes)
                # If we need to resample, go back somewhere into the past
                if bad_segments:
                    start_node = random.choice(bad_segments)
                    sm.elem_defs[start_node] = self.stat_source.sample_for(sm.bg, start_node)
                    built_nodes = built_nodes[:built_nodes.index(start_node)]
                    log.debug("Going back to node {}".format(start_node))
                else:
                    start_node = built_nodes[-1]
                    log.debug("Proceeding with next node {}".format(start_node))
                newbuilt_nodes = sm.new_traverse_and_build(start = start_node, max_steps = 1)
        except KeyboardInterrupt:
            print("No valid structure could be built after {} iterations".format(iterations), file=sys.stderr)
            raise
        for elem in _determined_broken_ml_segments(built_nodes, sm.bg):
            if elem in sm.junction_constraint_energy and hasattr(sm.junction_constraint_energy[elem], "used_stat"):
                used_stat = sm.junction_constraint_energy[elem].used_stat
                log.debug("Assigning stat %s to broken ml segment %s",
                          used_stat, elem)
                sm.elem_defs[elem] = used_stat
                sm.save_sampled_elems()
        log.debug("++++++++++++++++++++++++++++++++++++++")

    def _check_sampled_ml(self, sm, ml):
        """
        Raises an error if the sampled multiloop segment does not fulfill the junction closure energy.
        For the RoughJunctionClosureEnergy, this should only raise
        if the stat_source contains faulty stats.
        :param sm: The spatial model
        :param ml: The name of the junction segment we need to check (e.g. "m0")
        :raises: ValueError, if the ml-segment does not fulfill the energy
        :returns: None
        """
        if ml not in sm.junction_constraint_energy:
            return
        if sm.junction_constraint_energy[ml].eval_energy(sm.bg, nodes = [ml])!=0:
            dist = ftug.junction_virtual_atom_distance(sm.bg, ml)
            raise ValueError("Multiloop {} does not fulfill the constraints. "
                             "Sampled as {}, "
                             "distance = {}".format(ml, sm.elem_defs[ml], dist))

    def _get_bad_ml_segments(self, sm, nodes):
        """
        Return a list of bulges that are part of nodes and belog to a multiloop that doesn't
        fulfill the energy

        :param sm: The spatial model
        :param nodes: The built nodes of this spatial model. Only take these nodes and
                      fully defined ml-segments into account.
        :returns: A list of ml-elements that are part of a bad loop,
                  or an empty list, if the junction constriant energy is zero.
        """
        det_br_nodes = _determined_broken_ml_segments(nodes, sm.bg)
        log.debug("Evaluationg junction energy for nodes %s", det_br_nodes)
        for ml in det_br_nodes:
            if ml in sm.junction_constraint_energy:
                ej = sm.junction_constraint_energy[ml].eval_energy( sm.bg, nodes=det_br_nodes)
                log.debug("Junction Energy for nodes {} (=> {}) is {}".format(nodes, det_br_nodes, ej))
                if ej>0:
                    bad_loop_nodes =  [ x for x in nodes if x[0]=="m"]
                    log.debug("Bad loop nodes = {} (filtered from {})".format(bad_loop_nodes,
                                                                              sm.junction_constraint_energy[ml].bad_bulges))
                    return bad_loop_nodes
        return []

    def _get_bad_clash_segments(self, sm, nodes):
        """
        Return a list of interior loops and multiloop segments between the
        first stem in nodes that has a clash and the end of the structure.

        :param sm: The spatial model
        :param nodes: Only take these nodes into account
        :returns: A list of i and m element that were built after the first stem
                  with clashes, or an empty list is no clashes are detected.
        """
        if sm.constraint_energy is None:
            return []
        ec = sm.constraint_energy.eval_energy(sm.bg, nodes=nodes)
        log.debug("Clash Energy for nodes {} is {}".format(nodes, ec))
        if ec>0:
            bad_stems=set(x for clash_pair in sm.constraint_energy.bad_bulges for x in clash_pair)
            first = min(nodes.index(st) for st in bad_stems)
            assert first>=0
            clash_nodes =  [ x for x in nodes[first:] if x[0] in ["m", "i"]]
            log.debug("Clash nodes {}".format(clash_nodes))
            return clash_nodes
        return []

    def _rebuild_clash_only(self, sm, nodes, changable):
        """
        Tries to rebuild part of the structure to remove clashes.

        .. note::

            This is more efficient than self._build_with_energies if only clashes should be
            removed from a substructure because it avoids unnecessary energy evaluations.

        :param sm: The spatial model to build
        :param nodes: Take only these nodes into account for energy calculation
        :param changable: Only try to change on of these nodes. A list!
        :param tries: maximal tries before giving up and returning False

        :returns: True, if a clash_free structure was built.
        """
        if sm.constraint_energy is None:
            return True
        if not changable:
            return False
        for i in range(self.clash_only_tries):
            node = random.choice(changable)
            sm.elem_defs[node] = self.stat_source.sample_for(sm.bg, node)
            sm.new_traverse_and_build(start=node, end=nodes[-1])
            ec = sm.constraint_energy.eval_energy(sm.bg, nodes=nodes)
            if ec == 0:
                log.debug("_rebuild_clash_only for {} was successful after {} tries".format(nodes, i))
                return True
        log.debug("_rebuild_clash_only for {} was not successful.".format(nodes, i))
        return False

class FairBuilder(Builder):
    @profile
    def __init__(self, stat_source, output_dir = None, store_failed=False):
        """
        :param store_failed: Should structures, that do not fulfill the constraint energy, be stored?
                             A boolean or one of the following strings: "junction", "clash", "list"
                             In case of list: only append the failure reasons to the file clashlist.txt
        """
        super(FairBuilder, self).__init__(stat_source)
        self.output_dir = output_dir
        self.store_failed = store_failed
        self._failed_save_counter = 0
        self._success_save_counter = 0
    def build(self, sm):
        while True:
            self._attempt_to_build(sm)
            if self._fulfills_junction_energy(sm) and self._fulfills_clash_energy(sm):
                return

    def _attempt_to_build(self, sm):
        sm.sample_stats(self.stat_source)
        sm.new_traverse_and_build()

    def _fulfills_junction_energy(self, sm):
        if sm.fulfills_junction_energy():
            return True
        else:
            if self.store_failed is True or self.store_failed == "junction":
                self._store_failed(sm)
            elif self.store_failed=="list":
                with open(os.path.join(self.output_dir, "clashlist.txt"), "a") as f:
                    self._failed_save_counter += 1
                    bad_junctions = []
                    add = True
                    for j in self.junction_energy.bad_bulges:
                        if add:
                            bad_junctions.append(j)
                            add = False
                        else:
                            if j == bad_junctions[-1]:
                                add=True
                    f.write("{}: junction {}\n".format(self._failed_save_counter,
                                                     list(set(bad_junctions))))
            return False

    def _fulfills_clash_energy(self, sm):
        if sm.fulfills_clash_energy():
            return True
        else:
            if self.store_failed is True or self.store_failed == "clash":
                self._store_failed(sm)
            elif self.store_failed=="list":
                with open(os.path.join(self.output_dir, "clashlist.txt"), "a") as f:
                    self._failed_save_counter += 1
                    clash_pairs = set()
                    f.write("{}: clash {}\n".format(self._failed_save_counter, self.clash_energy.bad_bulges))
            return False

    def _store_failed(self, sm):
        self._failed_save_counter += 1
        with open(os.path.join(self.output_dir,
                              'failed{:06d}.coord'.format(self._failed_save_counter)), "w") as f:
            f.write(sm.bg.to_cg_string())

    def _store_success(self, sm):
        self._success_save_counter += 1
        with open(os.path.join(self.output_dir,
                              'build{:06d}.coord'.format(self._success_save_counter)), "w") as f:
            f.write(sm.bg.to_cg_string())

    @profile
    def success_probability(self, sm, target_attempts=None, target_structures=None, store_success = True):
        if target_attempts is None and target_structures is None:
            raise ValueError("Need target_structures or target_attempts")
        attempts = 0
        junction_failures = 0
        clashes = 0
        success = 0
        while True:
            if target_attempts is not None and attempts>=target_attempts:
                break
            self._attempt_to_build(sm)
            attempts+=1
            if not self._fulfills_junction_energy(sm):
                junction_failures += 1
                continue
            if not self._fulfills_clash_energy(sm):
                clashes += 1
                continue
            if store_success:
                self._store_success(sm)
            success+=1
            if target_structures is not None and success>target_structures:
                break
        log.info("Success_probability for fair building: {} attempts, thereof {} with failed "
                 "junctions, {} of the remaining structures have clashes, "
                 "{} were successful.".format(attempts, junction_failures, clashes, success))
        log.info("{:.0%} % junction failures, {:.0%}% clash, {:.0%}% "
                 " are ok.".format(junction_failures/attempts, clashes/attempts, success/attempts))
        good_j = attempts - junction_failures
        if good_j>0:
            log.info("For structures with good junctions: {:.0%}% clash, {:.0%}% "
                 " are ok.".format(clashes/good_j, success/good_j))
        return success, attempts, junction_failures, clashes

class ChangingMSTBuilder(FairBuilder):
    def _attempt_to_build(self, sm):
        if sm.bg.mst is None:
            sm.bg.traverse_graph()
        elem = random.choice(list(sm.bg.mst-sm.frozen_elements))
        if elem[0]=="m":
            try:
                sm.set_multiloop_break_segment(elem)
            except ValueError:
                pass
        sm.sample_stats(self.stat_source)
        sm.new_traverse_and_build()

class DimerizationBuilder(FairBuilder):
    #Inspired by dimerization for SAWs, as summerized in the review doi:10.1016/0920-5632(96)00042-4
    def __init__(self, stat_source, output_dir = None, store_failed=False):
        """
        Build all multiloops independently and then attempt to connect them
        via building of the rest of the structure. Reject the whole structure,
        if connecting of multiloops fails.

        :param store_failed: BOOL. Should structures, that do not fulfill the clash energy, be stored?
        :param output_dir: Where failed structures will be saved
        """
        super(DimerizationBuilder, self).__init__(stat_source, output_dir, store_failed)

    def _attempt_to_build(self, sm):
        sm.sample_stats(self.stat_source)
        sm.new_traverse_and_build()
        for multi_loop in sm.bg.find_mlonly_multiloops():
            log.debug("Sampling multiloop %r", multi_loop)
            while not self._fulfills_junction_energy(sm, multi_loop):
                self._change_multi_loop(sm, multi_loop)
        assert self._fulfills_junction_energy(sm)

    def _fulfills_junction_energy(self, sm, nodes=None):
        """
        In contrast to the super-class, this does not store failures,
        because we are only looking at partial structures here.
        """
        if nodes[0] in sm.junction_constraint_energy:
            return sm.junction_constraint_energy[nodes[0]]==0
        return True

    def _change_multi_loop(self, sm, multi_loop):
        node = random.choice(list(set(multi_loop) & sm.bg.mst))
        sm.elem_defs[node] = self.stat_source.sample_for(sm.bg, node)
        sm.new_traverse_and_build(start=node)


###############################################################################
###  Commandine arguments
###############################################################################

def update_parser(parser):
    builder_options = parser.add_argument_group("Building Spatial models",
                                    description="Options for specifying how to "
                                                "build the initial 3D structure.")
    builder_options.add_argument('-f', '--fair-building', action="store_true",
                                 help = "Try to build the structure using a fair \n"
                                        "but slow algorithm.\n "
                                        "This flag implies --start-from-scratch")
    builder_options.add_argument('--start-from-scratch', default=False, action='store_true',
                        help="Do not attempt to start at the input conformation.\n"
                             "(Automatically True for fasta files.)")

def from_args(args, stat_source):
    if args.fair_building:
        build_function = FairBuilder(stat_source).build
    else:
        b = Builder(stat_source)
        if args.start_from_scratch:
            build_function = b.build
        else:
            build_function = b.accept_or_build
    return build_function
