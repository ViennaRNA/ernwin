#!/usr/bin/python
from __future__ import print_function


from __future__ import absolute_import
import itertools as it
import random
import os.path as op
import numpy as np
import math
import sys
from collections import defaultdict
import warnings
import glob
import copy
import logging
import textwrap

import Bio.PDB as bpdb
import Bio.PDB.Chain as bpdbc

from logging_exceptions import log_to_exception
from commandline_parsable import split_by_outerlevel_character

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.stats as ftms
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.pdb as ftup
import forgi.threedee.utilities.vector as ftuv
import forgi.utilities.debug as fud

import fess.builder.config as cbc
import fess.builder.energy as fbe # Used for commandline-parsing
from fess.builder._commandline_helper import replica_substring
from six.moves import range


log = logging.getLogger(__name__)


class StemModel:
    '''
    A way of encapsulating the coarse grain 3D stem.
    '''

    def __init__(self, name=None, mids=None, twists=None):
        self.name = name

        if mids == None:
            self.mids = (np.array([0., 0., 0.]), np.array([0., 0., 1.0]))
        else:
            self.mids = mids
        if twists == None:
            self.twists = (np.array([0., 1., 0.]), np.array([1., 0., 0.0]))
        else:
            self.twists = twists

    def __str__(self):
        return str(self.mids) + '\n' + str(self.twists)

    def __eq__(self, other):
        mids0_close = np.allclose(self.mids[0], other.mids[0], atol=0.1)
        mids1_close = np.allclose(self.mids[1], other.mids[1], atol=0.1)
        twists0_close = np.allclose(self.twists[0], other.twists[0], atol=0.1)
        twists1_close = np.allclose(self.twists[1], other.twists[1], atol=0.1)

        return mids0_close and mids1_close and twists0_close and twists1_close

    def reverse(self):
        '''
        Reverse this stem's orientation so that the order of the mids
        is backwards. I.e. mids[1] = mids[0]...
        '''
        return StemModel(self.name, (self.mids[1], self.mids[0]), (self.twists[1], self.twists[0]))

    def vec(self, from_side_to_side = (0, 1)):
        from_side, to_side = from_side_to_side
        return self.mids[to_side] - self.mids[from_side]

    def rotate(self, rot_mat, offset=np.array([0., 0., 0.])):
        '''
        Rotate the stem and its twists according to the definition
        of rot_mat.
    assert np.array_equal(bm.mids[0], self.bg.coords[connected][1])

        @param rot_mat: A rotation matrix.
        '''
        self.mids = (np.dot(rot_mat, self.mids[0] - offset) + offset, np.dot(rot_mat, self.mids[1] - offset) + offset)
        self.twists = (np.dot(rot_mat, self.twists[0]), np.dot(rot_mat, self.twists[1]))

    def translate(self, translation):
        '''
        Translate the stem.
        '''
        self.mids = (self.mids[0] + translation, self.mids[1] + translation)

    def length(self):
        '''
        Get the length of this stem.
        '''
        return ftuv.magnitude(self.mids[1] - self.mids[0])

class BulgeModel:
    '''
    A way of encapsulating a coarse grain 3D loop.
    '''

    def __init__(self, mids=None):
        if mids == None:
            self.mids = (np.array([0., 0., 0.]), np.array([0., 0., 1.0]))
        else:
            self.mids = mids

    def __str__(self):
        return str(self.mids)

def translate_chain(chains, translation):
    '''
    Translate all of the atoms in a chain by a certain amount.

    :param chains: A dict chain_id:Bio.PDB.Chain instance to be translated.
    :translation: A vector indicating the direction of the translation.
    '''

    for chain in chains.values():
        atoms = bpdb.Selection.unfold_entities(chain, 'A')
        for atom in atoms:
            atom.transform(ftuv.identity_matrix, translation)

def rotate_chain(chains, rot_mat, offset):
    '''
    Move according to rot_mat for the position of offset.

    :param chains: A dict chain_id:Bio.PDB.Chain instance.
    :param rot_mat: A left_multiplying rotation_matrix.
    :param offset: The position from which to do the rotation.
    '''
    new_coords = []
    for chain in chains.values():

        atoms = bpdb.Selection.unfold_entities(chain, 'A')

        for atom in atoms:
            #atom.transform(ftuv.identity_matrix, -offset)
            assert ftuv.magnitude(atom.coord - offset) < ftuv.magnitude(atom.coord)
            atom.coord -= offset
            new_coords.append(atom.coord)
            atom.transform(rot_mat, offset)
        dev_from_cent = ftuv.magnitude(np.sum(new_coords, axis=0)/len(new_coords))
        if dev_from_cent>1:
            log.warning("{} not close to zero".format(dev_from_cent))



def get_stem_rotation_matrix(stem, stem2, use_average_method=False):
    """
    :param stem: The first StemModel
    :param stem2: The second StemModel

    :retuirns: A RotationMatrix.
               Use stem1.vec()*rotMat to rotate stem1 onto stem2
               Use rotMat*stem2.vec() to rotate stem2 onto stem1
    """
    #twist1 = (stem.twists[0] + stem.twists[1]) / 2.

    if not use_average_method:
        twist1 = stem.twists[0]
        twist2 = stem2.twists[0]

    else:
        twist1 = ftug.virtual_res_3d_pos_core(stem.mids, stem.twists, 2, 4)[1]
        twist2 = ftug.virtual_res_3d_pos_core(stem2.mids, stem2.twists, 2, 4)[1]

    return ftuv.get_double_alignment_matrix((stem.vec(), twist1),(stem2.vec(), twist2))
    # get normalvector to stem and twist.
    comp1 = np.cross(stem.vec(), twist1)

    # rotate around the first stem by t degrees
    rot_mat1 = ftuv.rotation_matrix(stem.vec(), t)
    rot_mat2 = ftuv.rotation_matrix(twist1, u - math.pi/2)
    rot_mat3 = ftuv.rotation_matrix(comp1, v)

    rot_mat4 = np.dot(rot_mat3, np.dot(rot_mat2, rot_mat1))

    return rot_mat4





def extract_stem_from_chain(chain, stem_def):
    '''
    Create a Chain consisting of just the atoms from the stem def.

    @param chain: The chain containing the stem (and potentially other atoms)
    @param stem_def: The define of the stem to be extracted
    '''
    c = bpdbc.Chain(' ')

    for i in range(stem_def.bp_length + 1):
        c.add(chain[stem_def.define[0] + i])
        c.add(chain[stem_def.define[2] + i])

    return c

def place_new_stem(prev_stem, stem_params, bulge_params, s1b_s1e, stem_name=''):
    '''
    Place a new stem with a particular orientation with respect
    to the previous stem.

    @param prev_stem: The already existing stem
    @param stem_params: A description of the new stem (StemStat)
    @param bulge_params: The AngleStat that specifies how the new stem will
                         be oriented with respect to the first stem.
    @param (s1b, s1e): Which side of the first stem to place the second on.
    '''
    s1b, s1e = s1b_s1e
    stem = StemModel()

    transposed_stem1_basis = ftuv.create_orthonormal_basis(prev_stem.vec((s1b, s1e)), prev_stem.twists[s1e]).transpose()
    log.debug("Place new stem: transposed_stem1_basis: %s", transposed_stem1_basis)
    start_location = ftug.stem2_pos_from_stem1_1(transposed_stem1_basis, bulge_params.position_params())
    log.debug("Start location: %s", start_location)
    stem_orientation = ftug.stem2_orient_from_stem1_1(transposed_stem1_basis, [stem_params.phys_length] + list(bulge_params.orientation_params()))
    log.debug("Stem_orientation: %s", stem_orientation)
    twist1 = ftug.twist2_orient_from_stem1_1(transposed_stem1_basis, bulge_params.twist_params())
    log.debug("twist1: %s", twist1)

    #assert np.allclose(np.dot(stem_orientation, twist1), 0)

    mid1 = prev_stem.mids[s1e] + start_location
    mid2 = mid1 + stem_orientation



    stem.mids = (mid1, mid2)

    log.debug("stem_params.twist_angle: %s", stem_params.twist_angle)
    twist2 = ftug.twist2_from_twist1(stem_orientation, twist1, stem_params.twist_angle)
    log.debug("twist2: %s", twist2)
    stem.twists = (twist1, twist2)

    #assert np.allclose(np.dot(stem_orientation, twist2), 0)
    return stem

def create_empty_energy():
    log.debug("Creating empty Energy for junction_constraint_energy-fdefaultdict")
    return fbe.CombinedEnergy()

class SpatialModel:
    '''
    A way of building RNA structures given angle statistics as well
    as length statistics.
    '''

    def __init__(self, bg, frozen_elements = None):
        '''
        Initialize the structure.

        @param bg: The BulgeGraph containing information about the coarse grained
                   elements and how they are connected.


        '''
        self.stems = dict()
        self.bulges = dict()

        if frozen_elements is None:
            self.frozen_elements = set()
        else:
            self.frozen_elements = set(frozen_elements)

        self.chain = bpdb.Chain.Chain(' ')
        self.build_chain = False
        self.constraint_energy = None
        #: A dictionary {broken_ml_segment: energy}
        self.junction_constraint_energy = defaultdict(create_empty_energy)

        self.elem_defs=dict()

        self.bg = bg
        # We plan to modify the structure, so discard the cahin.
        # This avoids a bug with virtual residues of loops.
        self.bg.chains=None
        self.add_to_skip()

        try:
            self.bg.add_all_virtual_residues()
        except (ftmc.RnaMissing3dError, AssertionError):
            # The structure is probably new and doesnt have coordinates yet
            pass

    def sample_stats(self, stat_source):

        for d in self.bg.defines:
            #Frozen elements are only sampled once, if they couldn't be loaded.
            if d in self.frozen_elements and d in self.elem_defs:
                continue
            elif d in self.frozen_elements:
                log.warning("Frozen element %s sampled once, because no stats were present.", d)
            if d[0] == 'm':
                if self.bg.get_angle_type(d, allow_broken = False) is None:
                    # this section isn't sampled because a multiloop
                    # is broken here
                    continue
            try:
                self.elem_defs[d] = stat_source.sample_for(self.bg, d)
            except:
                log.error("Error sampling stats for element %s.", d)
                raise

    def sampled_from_bg(self):
        '''
        Get the information about the sampled elements from the underlying BulgeGraph.
        '''
        # get the stem defs
        # self.get_sampled_bulges()

        raise Exception("This needs to be re-written, possible using FilteredConformationStats")


    def create_native_stem_models(self):
        '''
        Create StemModels from the stem definitions in the graph file.
        '''
        stems = dict()

        for d in self.bg.defines.keys():
            if d[0] == 's':
                stems[d] = StemModel(d, self.bg.coords[d], self.bg.twists[d])

                if self.build_chain:
                    reconstruct_stem(self, d, self.chain, stem_library=cbc.Configuration.stem_library, stem=stems[d])

        self.stems = stems


    def add_loop(self, name, prev_stem_node, params=None, loop_defs=None):
        '''
        Connect a loop to the previous stem.
        '''
        if loop_defs == None:
            loop_defs = self.elem_defs

        prev_stem = self.stems[prev_stem_node]
        (s1b, s1e) = self.bg.get_sides(prev_stem_node, name)

        if params == None:
            r = loop_defs[name].phys_length
            u = loop_defs[name].u
            v = loop_defs[name].v

            params = (r, u, v)

        start_mid = prev_stem.mids[s1b]
        log.debug("prev stem %s: mids[%s]=%s, coords = %s", prev_stem_node, s1b, prev_stem.mids[s1b], self.bg.coords[prev_stem_node])
        (r, u, v) = params

        direction = ftug.stem2_pos_from_stem1(prev_stem.vec((s1e, s1b)), prev_stem.twists[s1b], (r, u, v))
        end_mid = start_mid + direction
        log.debug("loop added from %s to %s", start_mid, end_mid)
        self.bulges[name] = BulgeModel((start_mid, end_mid))

    def find_start_node(self):
        '''
        Find a node from which to begin building. This should ideally be a loop
        region as from there we can just randomly orient the first stem.
        '''

        edge = next(self.bg.sorted_stem_iterator())
        define = 'start'
        return (edge, define, StemModel(edge))


    def save_sampled_elems(self):
        '''
        Save the information about all of the sampled elements.
        '''
        for d, ed in self.elem_defs.items():
            try:
                self.bg.sampled[d] = [ed.pdb_name] + [len(ed.define)] + ed.define
                self.bg.vbase[d] = ed.vbase
                self.bg.vsugar[d] = ed.vsugar
                self.bg.vbackbone[d] = ed.vbackbone
                if d[0]!="s":
                    try:
                        log.debug("Updating vposs from %s. Were %s", d, self.bg.vposs[d] )
                    except:
                        pass
                    self.bg.vposs[d]=ed.vres
                    log.debug("Set %s to %s ", d, ed.vres)
                if d not in self.bg.get_mst():
                    self.bg.infos["vstat_{}".format(d)]=[str(ed)]
            except:
                log.debug("Error setting {}".format(d))
                raise

    def change_mst(self, new_mst, stat_source):
        """
        Set bg.mst to new_mst and recalculate invalid values.

        Assumes new_mst is a valid minimal spanning tree.
        """
        old_only = self.bg.mst - new_mst
        self.bg.mst = new_mst

        self.bg.build_order = None #No longer valid
        self.bg.ang_types = None
        self.load_sampled_elems(stat_source=None)
        #self.new_traverse_and_build(start='start', include_start = True)

    def set_multiloop_break_segment(self, d):
        """
        Change the minimum spanning tree of the Coarse grain RNA in a way that
        *  The multiloop segment d is broken.
        *  The mst is still connected.
        *  The coordinates of all coarse grain elements are not changed if `self.traverse_and_build` is called.

        :param d: The ML segment that should be broken. E.g. "m1"
        :returns: The ML segment that was broken before but now was added to the build_order.
        """
        log.info("Changing MST")
        if self.bg.mst is None:
            self.bg.traverse_graph()
        junction_nodes = set(self.bg.shortest_mlonly_multiloop(d))
        missing_nodes = junction_nodes - self.bg.mst
        if d in missing_nodes:
            return None #The specified cg element is already a breaking point.
        if len(missing_nodes)==1: #The easy case. Just exchange the two ml segments
            log.info("One ml-segment missing: %s", missing_nodes)
            log.info("mst: %s", self.bg.mst)
            self.bg.mst.remove(d)
            log.info("mst: %s", self.bg.mst)
            self.bg.mst|=missing_nodes
            log.info("mst: %s", self.bg.mst)

            self.bg.build_order = None #No longer valid
            self.bg.ang_types = None
            if self.elem_defs and d in self.elem_defs:
                del self.elem_defs[d]
            self.bg.traverse_graph()
            log.info("mst: %s", self.bg.mst)

            self.load_sampled_elems(stat_source=None)
            new_node, = missing_nodes
            return new_node
        elif len(missing_nodes)>=2:
            log.info("More than one ml-segment missing")

            #Delete requested Edge
            self.bg.mst.remove(d)

            #Find connected parts of the structure.
            forest = []
            print("MST", self.bg.mst)
            for m in self.bg.mst:
                neighbors = list(self.bg.edges[m])
                m_and_s = [m]
                for n in neighbors:
                    if n in self.bg.mst:
                        m_and_s.append(n)
                for tree in forest:
                    if any(loop in tree for loop in m_and_s):
                        tree |= set(m_and_s)
                        break
                else:
                    forest.append(set(m_and_s))
                while True: #I do not know if a loop+2 stems can ever connect more than 2 trees. Maybe not needed????
                    for t1, t2 in it.combinations(forest, 2):
                        if len(t1 & t2) > 0: #Non-empty intersection
                            print("t1&t2:", t1, t2, t1&t2)
                            new_tree = t1 | t2 #Union Tree
                            forest.remove(t1)
                            forest.remove(t2)
                            forest.append(new_tree)
                            break
                    else:
                        break
            assert len(forest)==2 # No matter how many missing nodes, this should always be 2
            for missing_node in missing_nodes:
                neighbors = list(self.bg.edges[missing_node])
                if set(neighbors) & forest[0] and set(neighbors) & forest[1]:
                    self.bg.mst.add(missing_node)
                    new_node=missing_node
                    break
            else:
                raise ValueError("Cannot break loop {:0}. Cannot connect RNA if {:0} is broken.".format(d))
            self.bg.build_order = None #No longer valid
            self.bg.ang_types = None
            if self.elem_defs and d in self.elem_defs:
                del self.elem_defs[d]
            self.bg.traverse_graph()
            self.load_sampled_elems(stat_source=None)
            return new_node
        elif len(missing_nodes)==0:
            raise ValueError("Cannot break loop {}. This multi loop is not cyclic and thus has no broken fragment.".format(d))

    def load_sampled_elems(self, stat_source):
        '''
        Load information from the CoarseGrainRNA into self.elem_defs

        ..note: This is called by the Builder to avoid building the structure from scratch.
        '''
        build_order = self.bg.traverse_graph()
        for d in self.bg.defines.keys():
            if d[0]=="s":
                stat=self.bg.get_stem_stats(d)
            elif d[0] in ("h","f","t"):
                stat=self.bg.get_loop_stat(d)
                #stat.define=self.bg.defines[d]
            elif d[0] in ("m", "i"):
                stat=ftms.AngleStat()
                #load angle_stats in direction of build order!
                for bo in build_order:
                    if bo[1]==d:
                        if self.bg.defines[bo[0]]<self.bg.defines[bo[2]]:
                            forward=True
                        else:
                            forward=False
                        stat=self.bg.get_bulge_angle_stats_core(d, forward=forward)
                        break
                else: #Not in build_order. Probably broken ml segment
                    assert d[0]=="m"
                    #If it was frozen or stored in the file, we add its stats to elem_defs.
                    if d in self.frozen_elements:
                        stat,=[s for s in self.bg.get_bulge_angle_stats(d) if s.ang_type==self.bg.get_angle_type(d, allow_broken=True)]
                        log.info("Loading stat from cg for frozen ml segment %s to %s", d, stat)
                    else:
                        #Do not add this element to elem_defs!
                        log.info("Not adding stat to broken ml segment %s from cg file", d)
                        continue
            log.debug("Loading stat %s from bg to elem_defs[%s]", stat, d)
            self.elem_defs[d]=stat
        if stat_source is not None:
            if self.bg.sampled:
                log.info("Now overwriting stats by sampled stats from stat_source")
            build_order = self.bg.traverse_graph()
            for elem, sampled in self.bg.sampled.items():
                pdb_name=sampled[0]
                stat = stat_source.load_stat_by_name(self.bg, elem, pdb_name)
                self.elem_defs[elem]=stat

    def get_transform(self, edge):
        '''
        Get the location of one end of a stem.

        Used in aligning a group of models around a single edge.

        @param edge: The name of the edge.
        '''
        assert(edge[0] == 's')

        return self.stems[edge].mids[0]

    def get_rotation(self, edge):
        '''
        Get a rotation matrix that will align a point in the direction of a particular
        edge.

        Used in aligning a group of models around a single edge.

        @param edge: The name of the target edge.
        '''
        assert(edge[0] == 's')

        target = [np.array([1., 0., 0.]), np.array([0., 1., 0.])]


        vec1 = self.stems[edge].vec()
        twist1 = self.stems[edge].twists[1]

        mat = ftuv.get_double_alignment_matrix(target, [vec1, twist1])

        return mat

    def add_stem(self, stem_name, stem_params, prev_stem, bulge_params, s1b_s1e):
        '''
        Add a stem after a bulge.

        The bulge parameters will determine where to place the stem in relation
        to the one before it (prev_stem). This one's length and twist are
        defined by the parameters stem_params.

        @param stem_params: The parameters describing the length and twist of the stem.
        @param prev_stem: The location of the previous stem
        @param bulge_params: The parameters of the bulge.
        @param side: The side of this stem that is away from the bulge
        '''
        s1b, s1e = s1b_s1e
        log.debug("add_stem {}".format(stem_name))
        stem = place_new_stem(prev_stem, stem_params, bulge_params, (s1b, s1e), stem_name)

        stem.name = stem_name

        if self.build_chain:
            reconstruct_stem(self, stem_name, self.chain, stem_library=cbc.Configuration.stem_library, stem=stem)

        return stem

    def fill_in_bulges_and_loops(self):
        log.debug("Started fill_in_bulges_and_loops")
        loops = list(self.bg.hloop_iterator())
        fiveprime = list(self.bg.floop_iterator())
        threeprime = list(self.bg.tloop_iterator())

        for d in self.bg.defines.keys():
            if d[0] != 's':
                if d in loops or d in fiveprime or d in threeprime:
                    log.debug("Adding loop {} (connected to {})".format(d, list(self.bg.edges[d])[0]))
                    self.add_loop(d, list(self.bg.edges[d])[0])
                else:
                    connections = self.bg.connections(d)

                    # Should be a bulge connecting two stems
                    assert(len(connections) == 2)
                    for conn in connections:
                        assert(conn[0] == 's')

                    (s1b, s1e) = self.bg.get_sides(connections[0], d)
                    (s2b, s2e) = self.bg.get_sides(connections[1], d)

                    s1mid = self.stems[connections[0]].mids[s1b]
                    s2mid = self.stems[connections[1]].mids[s2b]

                    self.bulges[d] = BulgeModel((s1mid, s2mid))

    def stem_to_coords(self, stem):
        """

        """
        log = logging.getLogger(__name__+".SpatialModel.stem_to_coords")
        sm = self.stems[stem]
        if log.isEnabledFor(logging.DEBUG):
            if (not np.allclose(self.bg.coords[stem][0], sm.mids[0]) or
                not np.allclose(self.bg.coords[stem][1], sm.mids[1])):
                log.debug("Changing stem %s: %s to %s",stem, self.bg.coords[stem], (sm.mids[0], sm.mids[1]))
            if (not np.allclose(self.bg.twists[stem][0], sm.twists[0]) or
                not np.allclose(self.bg.twists[stem][1], sm.twists[1])):
                log.debug("Changing stem twist %s : %s to %s", stem, self.bg.twists[stem], (sm.twists[0], sm.twists[1]))

        self.bg.coords[stem] = np.array([sm.mids[0], sm.mids[1]])
        self.bg.twists[stem] = np.array([sm.twists[0], sm.twists[1]])

    def _loops_to_coords(self):
        '''
        Add all of the stem and bulge coordinates to the BulgeGraph data structure.
        '''

        log.debug("_loops_to_coords: Adding bulge coodinates from stems")
        self.bg.add_bulge_coords_from_stems()

        for d in self.bg.hloop_iterator():
            bm = self.bulges[d]
            connected, =self.bg.edges[d]
            if not np.allclose(bm.mids[0], self.bg.coords[connected][1]):
                log.error("BulgeModel for {} with coords {} and {}".format(d, bm.mids[0], bm.mids[1]))
                log.error("Connected to stem {} with coords {} and {}".format(connected, self.bg.coords[connected][0], self.bg.coords[connected][1]))
                assert False, "Bulge {}, Difference {}".format(d, bm.mids[0]-self.bg.coords[connected][1])
            self.bg.coords[d] = np.array([bm.mids[0], bm.mids[1]])
        for d in it.chain(self.bg.floop_iterator(), self.bg.tloop_iterator()):
            if d in self.bg.defines:
                bm = self.bulges[d]
                connected, =self.bg.edges[d]
                # We do not know at which stem-side the 3'/5' element is, because of cofold structures.
                assert np.array_equal(bm.mids[0], self.bg.coords[connected][0]) or np.array_equal(bm.mids[0], self.bg.coords[connected][1])
                self.bg.coords[d] = np.array([bm.mids[0], bm.mids[1]])

    def get_sampled_bulges(self):
        '''
        Do a breadth first traversal and return the bulges which are
        sampled. This will be used to determine which ones are closed.
        '''
        visited = set()
        prev_visited = set()

        first_node = self.find_start_node()[:2]
        self.sampled_bulges = []
        self.visit_order = []

        to_visit = [first_node]

        while len(to_visit) > 0:
            to_visit.sort(key=lambda x: -self.bg.stem_length(x[0]))
            #rand.shuffle(to_visit)
            (curr_node, prev_node) = to_visit.pop()

            while curr_node in visited:
                if len(to_visit) > 0:
                    (curr_node, prev_node) = to_visit.pop()
                else:
                    self.visit_order = visited
                    return self.sampled_bulges

            #print curr_node, prev_node

            visited.add(curr_node)
            prev_visited.add(prev_node)

            if curr_node[0] == 's':
                self.sampled_bulges += [prev_node]

            for edge in self.bg.edges[curr_node]:
                if edge not in visited:
                    to_visit.append((edge, curr_node))

        self.visit_order = visited
        return self.sampled_bulges

        #self.prev_visit_order = prev_visited

    def _finish_building(self):
        log.debug("Finish building")
        log.debug("(1) vposs now %s", self.bg.vposs)
        self.fill_in_bulges_and_loops()
        log.debug("(2) vposs now %s", self.bg.vposs)
        self._loops_to_coords()
        log.debug("(3) vposs now %s", self.bg.vposs)
        self.save_sampled_elems()
        log.debug("(4) vposs now %s", self.bg.vposs)
        self.bg.add_all_virtual_residues()
        log.debug("(5) vposs now %s", self.bg.vposs)

    def add_to_skip(self):
        '''
        Build a minimum spanning tree of the bulge graph.
        '''
        to_skip = set()
        visited = set()

        for d in self.bg.defines.keys():
            if d in visited:
                continue

            visited.add(d)

            loop = self.bg.find_bulge_loop(d, 1000000) + [d]
            if len(loop) > 1:
                loop_w_sizes = [(self.bg.stem_length(l), l) for l in loop if l[0] != 's']
                loop_w_sizes += [(0, l) for l in loop if l[0] == 's']
                to_remove = max(loop_w_sizes)[1]
                to_skip.add(to_remove)

            for l in loop:
                visited.add(l)

        self.to_skip = to_skip


    def new_traverse_and_build(self, start='start', max_steps=float('inf'),
                               end=None, include_start=False, finish_building = True):
        '''
        Build a 3D structure from the graph in self.bg and the stats from self.elem_defs.
        :param start: Optional; Start building the given element. If it is a stem, build AFTER this stem.
                      TODO: instead of linear build_order use a tree. Only build leaves dependend on the changed node

        :param max_staps: Optional; Build at most that many stems.
        :param end: Optional; End building once the given node is built.
                    If `end` and `max_steps` are given, the criterion that kicks in earlier counts.
        :param include_start: If True, build including the node given as start,
                              if False, only build AFTER it. Only has an effect is start is a stem.
        :param finish_building: Usually True. Set it to False, if self._finish_building
                                should not be called after the stems have been built.
        :returns: A list of course_grained elements that have been built.
                  (Useful, if start or max_steps is given).
        '''
        log.debug("new_traverse_and_build(self, start={}, max_steps={}, end={})".format(start, max_steps, end))
        if not self.bg.build_order:
            build_order = self.bg.traverse_graph()
        else:
            build_order = self.bg.build_order
        log.debug("build_order: %s", build_order)
        def buildorder_of(stemid, include = False):
            """
            Returns the buildorder of the multi-/ interior loop before the stem with stemid.
            :param stemid: a string describing a stem or loop, e.g. 's0', 'i3'
            :param include: Whether or not to include the starting node (only for stems),
                            for bulges, the stem AFTER the bulge is always the first stem placed.
            """

            if stemid == "s0":
                return int(include)

            if stemid.startswith('s'):
                for i, stem_loop_stem in enumerate(build_order):
                    if stemid==stem_loop_stem[2]:
                        if include:
                            return i
                        else:
                            return i+1
            else:
                if stemid[0] in "fth":
                    return float("inf") #Only fill in bulges and loops if loop was changed
                for i, stem_loop_stem in enumerate(build_order):
                    if stemid==stem_loop_stem[1]:
                        return i
            raise ValueError("{} not found in {}.".format(stemid,build_order))

        nodes = []
        if start == "start" or (start == "s0" and include_start):
            # add the first stem in relation to a non-existent stem
            first_stem = "s0"
            log.debug("new_traverse_and_build: Setting self.stems[{}] (=first  stem)".format(first_stem))
            self.stems[first_stem] = self.add_stem(first_stem, self.elem_defs[first_stem], StemModel(),
                                                   ftms.AngleStat(), (0,1))
            self.stem_to_coords(first_stem)
            nodes.append(first_stem)
            build_step = 0
        elif start=="end" or start[0] in "fth":
            self._finish_building()
            return []
        else:
            build_step = buildorder_of(start, include_start)
            if build_step >= len(build_order):
                if finish_building:
                    self._finish_building()
                return []
            prev_stem = build_order[build_step][0]
            try:
                log.debug("new_traverse_and_build: Checking self.stems[{}] (=prev_stem)".format(prev_stem))
                self.stems[prev_stem]
            except KeyError:
                raise ValueError("Cannot build structure starting from {0}, because the parts "
                                 "of the structure before {0} have never been built. "
                                 "(The start-option is only for RE-building)".format(start))

        max_build_steps = min(build_step+max_steps, len(build_order))
        while build_step < max_build_steps:
            (s1, l, s2) = build_order[build_step]
            nodes += [l, s2]
            build_step +=1

            prev_stem = self.stems[s1]
            angle_params = self.elem_defs[l]
            stem_params = self.elem_defs[s2]
            ang_type = self.bg.connection_type(l, [s1,s2])
            connection_ends = self.bg.connection_ends(ang_type)

            # get the direction of the first stem (which is used as a
            # coordinate system)
            if connection_ends[0] == 0:
                (s1b, s1e) = (1, 0)
            elif connection_ends[0] == 1:
                (s1b, s1e) = (0, 1)

            log.debug("new_traverse_and_build: Setting self.stems[{}] (connected to {} via {})".format(s2, s1, l))
            #log.debug("angle_params {}, stem_params {}, ang_type {}, connection_ends {}".format(angle_params, stem_params, ang_type, connection_ends))
            #log.debug("prev. stem MIDS: {}, TWISTS: {}".format(prev_stem.mids, prev_stem.twists))

            stem = self.add_stem(s2, stem_params, prev_stem,
                                 angle_params, (s1b, s1e))

            # check which way the newly connected stem was added
            # if its 1-end was added, the its coordinates need to
            # be reversed to reflect the fact it was added backwards
            if connection_ends[1] == 1:
                self.stems[s2] = stem.reverse()
            else:
                self.stems[s2] = stem

            self.stem_to_coords(s2)

            #Optional end-criterion given as a node label.
            if end is not None and end in nodes:
                break
        if finish_building:
            self._finish_building()
        return nodes


    def ml_stat_deviation(self, ml, stat):
        """
        Calculate the deviation in angstrom between the stem that would be placed using the given
        stats for an open multiloop segment and the true stem position.

        :param ml: A element name, e.g. "m0". Should correspond to a broken multiloop.
        :param stat: The fictive stat to use for this ml segment.
        """
        stat1, stat2 = self.bg.get_stats(ml)

        if stat.ang_type == stat1.ang_type:
            virtual_stat = stat1
        else:
            assert stat.ang_type == stat2.ang_type
            virtual_stat = stat2
        diff = stat.diff(virtual_stat, next_stem_length = 3)
        return diff

    def fulfills_constraint_energy(self):
        return self.fulfills_clash_energy() and self.fulfills_junction_energy()
    def fulfills_junction_energy(self):
        for mloop in self.bg.find_mlonly_multiloops():
            for loop in mloop:
                log.debug("Trying junction constraint energy for %s. Energies are %s", mloop, list(self.junction_constraint_energy.keys()))
                if loop not in self.bg.get_mst() and loop in self.junction_constraint_energy:
                    log.debug("Evaluating junction constraint energy %s", self.junction_constraint_energy[loop].shortname)
                    if self.junction_constraint_energy[loop].eval_energy(self.bg, nodes=mloop, sampled_stats=self.elem_defs)>0:
                        log.info("Junction {} is not closed".format(mloop))
                        return False
        return True
    def fulfills_clash_energy(self):
        if self.constraint_energy is None:
            warnings.warn("Model has no clash energy!")
        if self.constraint_energy is not None and self.constraint_energy.eval_energy(self.bg)>0:
            log.info("CLASHING")
            return False
        return True

#############################################################################################
### COMMANDLINE OPTIONS
#############################################################################################
def _auto_contrib_for_loop(loop, sm, stat_source):
    """
    Choose the energy based on the number of available stats.

    For non-regular multiloops or if only few stats are available for the
    broken ml-segments or if the total number of stat-combinations for
    the whole loop is too little, we use JDIST, otherwise we use M8[1FJC1]
    """
    if "regular_multiloop" not in sm.bg.describe_multiloop(loop):
        log.info("AUTO-contrib for non-regular multiloop %s is JDIST", loop)
        return "{}:JDIST".format(loop[0])
    broken_mls = set(loop)-sm.bg.get_mst()
    for ml in broken_mls:
        num_stats = sum(1 for _ in stat_source.iterate_stats_for(sm.bg, ml))
        if num_stats<50:
            log.info("AUTO-contrib for multiloop %s where "
                     "elem %s has only %s stats is JDIST", loop, ml, num_stats)
            return "{}:JDIST".format(ml)
    prod = 1
    for ml in loop:
        num_stats = sum(1 for _ in stat_source.iterate_stats_for(sm.bg, ml))
        prod *= num_stats
    if prod<10**6:
        log.info("AUTO-contrib for multiloop %s with "
                 "%e stat combinations is JDIST", loop, prod)
        return "{}:JDIST".format(ml)
    else:
        log.info("AUTO-contrib for multiloop %s is MAX8[1FJC1]", loop)
        return "{}:MAX8[1FJC1]".format(ml)

def _perml_energy_to_sm(sm, energy_string, stat_source):
    log.info("Setting up constraint energy from string %s", energy_string)
    all_loops = sm.bg.find_mlonly_multiloops()
    if not all_loops:
        log.info("Ignoring perml-energy. No loops present")
        return
    if energy_string == "AUTO":
        contribs = []
        for loop in  all_loops:
            contribs.append(_auto_contrib_for_loop(loop, sm, stat_source))
        energy_string = ",".join(contribs)
    mst=sm.bg.get_mst()
    if energy_string and energy_string != 'N':
        contributions = split_by_outerlevel_character(energy_string, ",")
        for contrib in contributions:
            parts = split_by_outerlevel_character(contrib, ":")
            if len(parts)==2:
                elem, energy_string = parts
            elif len(parts)==1:
                elem = None
                energy_string = contrib
            else:
                raise ValueError("Too many colons in energy specification `{}` "
                                 "(at most 1 allowed)".format(contribution))
            energy, = fbe.EnergyFunction.from_string(energy_string, cg=sm.bg,
                                                     iterations=None,
                                                     stat_source=stat_source,
                                                     elements=elem)
            if not hasattr(energy, "can_constrain") or energy.can_constrain != "junction":
                e = ValueError("The energy '{}' cannot be used as a junction "
                                 "constraint energy".format(energy_string))
                with log_to_exception(log, e):
                    if hasattr(energy, "can_constrain"):
                        log.error("can_constrain is %r", energy.can_constrain)
                    else:
                        log.error("Energy of type %r has no attr 'can_constrain'", type(energy).__name__)
                raise e
            log.info("Searching for loops that contain element %s", elem)
            for loop in all_loops:
                log.info("Potentially assigning something to %s", loop)
                if not elem or elem in loop:
                    for loop_elem in loop:
                        if loop_elem in mst:
                            continue
                        if hasattr(energy, "element"):
                            energy, = fbe.EnergyFunction.from_string(energy_string, cg=sm.bg, iterations=None, stat_source=stat_source, element=loop_elem)
                        log.info("Assigning %s to %s", energy.shortname, loop_elem)
                        sm.junction_constraint_energy[loop_elem].energies.append(energy)

def update_parser(parser):
    sm_options = parser.add_argument_group("Options related to the Spatial Model",
                                    description="These options control the SpatialModel")
    sm_options.add_argument('--freeze', type=str, default="",
                            help= "A comma-seperated list of cg-element names.\n"
                                  "These elements will not be changed during samplig.")
    sm_options.add_argument('--mst-breakpoints', type=str,
                            help="During initial MST creation, prefer to \n"
                                 "break the multiloops at the indicated nodes.\n"
                                 "A comma-seperated list. E.g. 'm0,m10,m12'")
    sm_options.add_argument('--constraint-energy-clash', type=str, default="CLASH",
                            help="Specify the constraint energies that require the complete\n"
                                 "spatial model for evaluation. Example: The clash energy. \n"
                                 "Use N for no energy")
    sm_options.add_argument('--constraint-energy-per-ml', type=str, default="AUTO",
                            help=textwrap.dedent("""\
                                    Specify the constraint energies that can be evaluated on the
                                    junction level without need for the complete spatial model.
                                    `AUTO` to automatically choose between JDIST and M?SFJ
                                    depending on the junction characteristics.
                                    `N` or empty string for no energy.
                                    Specify energies as documented for the `--energy` option.
                                    Normally, these energies apply to all junctions.
                                    Prefix contributions by a multiloop segment name followed by
                                    a colon, to apply the following contribution only to the
                                    multiloop containing this segment.
                                    Example: `m1:JDIST,M8[1FJC1]` applies the JDIST energy only
                                    to the junction containing the m1 segment, and additionally
                                    the M8[1SFJ1] energy to all multiloops.
                                    """))

def some_replica_different(args):
    """
    Return True, if the user specified
    different energies for different replicas.
    """
    return ('@' in args.constraint_energy_clash or
            '@' in args.constraint_energy_per_ml)

def from_args(args, cg, stat_source, replica=None):
    frozen = set(args.freeze.split(","))
    sm = SpatialModel(cg, frozen_elements=frozen )

    if args.mst_breakpoints:
        bps = args.mst_breakpoints.split(",")
        for bp in bps:
            sm.set_multiloop_break_segment(bp)

    if args.constraint_energy_clash != "N":
        clash_string = replica_substring(args.constraint_energy_clash, replica)
        clash_e = fbe.EnergyFunction.from_string( clash_string, cg=cg, iterations=None)
        sm.constraint_energy = fbe.CombinedEnergy(clash_e)
    if not sm.constraint_energy:
        log.error("WARNING: Not using constraint energy for SM")
    perml_string = replica_substring(args.constraint_energy_per_ml, replica)
    _perml_energy_to_sm(sm, perml_string, stat_source)

    return sm
