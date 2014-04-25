#!/usr/bin/python

import Bio.PDB as bpdb
import Bio.PDB.Chain as bpdbc
import itertools as it
import random
import os
import os.path as op
import warnings
import numpy as np
import numpy.linalg as nl
import math
import sys
import collections as c
import random as rand
import collections as c

import fess.builder.config as cbc
import forgi.threedee.model.stats as cbs
import forgi.threedee.utilities.graph_pdb as cgg
import forgi.threedee.utilities.pdb as ftup
import forgi.threedee.utilities.vector as cuv
import forgi.utilities.debug as fud
import forgi.threedee.model.coarse_grain as ftmc

from random import choice, uniform
from math import pi

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

    def vec(self, (from_side, to_side) = (0, 1)):
        return self.mids[to_side] - self.mids[from_side]

    def rotate(self, rot_mat, offset=np.array([0., 0., 0.])):
        '''
        Rotate the stem and its twists according to the definition
        of rot_mat.

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
        return cuv.magnitude(self.mids[1] - self.mids[0])

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
            
def translate_chain(chain, translation):
    '''
    Translate all of the atoms in a chain by a certain amount.

    @param chain: A Bio.PDB.Chain instance to be translated.
    @translation: A vector indicating the direction of the translation.
    '''
    atoms = bpdb.Selection.unfold_entities(chain, 'A')

    for atom in atoms:
        atom.transform(cuv.identity_matrix, translation)

def rotate_chain(chain, rot_mat, offset):
    '''
    Move according to rot_mat for the position of offset.

    @param chain: A Bio.PDB.Chain instance.
    @param rot_mat: A left_multiplying rotation_matrix.
    @param offset: The position from which to do the rotation.
    '''

    atoms = bpdb.Selection.unfold_entities(chain, 'A')

    for atom in atoms:
        #atom.transform(np.eye(3,3), -offset)
        atom.coord -= offset
        atom.transform(rot_mat, offset)

def define_to_stem_model(cg, chain, define):
    '''
    Extract a StemModel from a Bio.PDB.Chain structure.

    The define is 4-tuple containing the start and end coordinates
    of the stem on each strand. 

    s1s s1e s2s s2e

    @param chain: The Bio.PDB.Chain representation of the chain
    @param define: The BulgeGraph define
    @return: A StemModel with the coordinates and orientation of the stem.
    '''
    stem = StemModel(name=define)

    mids = cgg.get_mids(cg, chain, define)
    #mids = cgg.estimate_mids_core(chain, int(define[0]), int(define[3]), int(define[1]), int(define[2])) 

    stem.mids = tuple([m.get_array() for m in mids])
    stem.twists = cgg.get_twists(cg, chain, define)

    return stem

def get_stem_rotation_matrix(stem, (u, v, t), use_average_method=False):
    #twist1 = (stem.twists[0] + stem.twists[1]) / 2.

    if not use_average_method:
        twist1 = stem.twists[0]
    else:
        twist1 = cgg.virtual_res_3d_pos_core(stem.mids, stem.twists, 2, 4)[1]

    # rotate around the stem axis to adjust the twist

    # rotate down from the twist axis
    comp1 = np.cross(stem.vec(), twist1)

    rot_mat1 = cuv.rotation_matrix(stem.vec(), t)
    rot_mat2 = cuv.rotation_matrix(twist1, u - math.pi/2)
    rot_mat3 = cuv.rotation_matrix(comp1, v)

    rot_mat4 = np.dot(rot_mat3, np.dot(rot_mat2, rot_mat1))

    return rot_mat4

def align_chain_to_stem(cg, chain, define, stem2, use_average_method=False):
    stem1 = define_to_stem_model(cg, chain, define)
    tw1 = cgg.virtual_res_3d_pos_core(stem1.mids, stem1.twists, 2, 4)[1]
    tw2 = cgg.virtual_res_3d_pos_core(stem2.mids, stem2.twists, 2, 4)[1]

    '''
    (r, u, v, t) = cgg.get_stem_orientation_parameters(stem1.vec(), 
                                                       (stem1.twists[0] + stem1.twists[1]) / 2., 
                                                       stem2.vec(), 
                                                       (stem2.twists[0] + stem2.twists[1]) / 2.)
    '''
    if not use_average_method:
        (r, u, v, t) = cgg.get_stem_orientation_parameters(stem1.vec(), 
                                                           stem1.twists[0], 
                                                           stem2.vec(), 
                                                           stem2.twists[0])
    else:
        (r, u, v, t) = cgg.get_stem_orientation_parameters(stem1.vec(), 
                                                           tw1, 
                                                           stem2.vec(), 
                                                           tw2)
    rot_mat = get_stem_rotation_matrix(stem1, (math.pi-u, -v, -t), use_average_method)
    rotate_chain(chain, np.linalg.inv(rot_mat), (stem1.mids[0] + stem1.mids[1]) / 2.)
    translate_chain(chain, (stem2.mids[0] + stem2.mids[1]) / 2. - (stem1.mids[0] + stem1.mids[1]) / 2.)

def reconstruct_stem_core(cg_orig, stem_def, orig_def, new_chain, stem_library=dict(), stem=None, use_average_method=True):
    '''
    Reconstruct a particular stem.
    '''
    pdb_filename = op.expanduser(op.join('~/doarse/', stem_def.pdb_name, "temp.pdb"))
    cg_filename = op.expanduser(op.join('~/doarse/', stem_def.pdb_name, "temp.cg"))

    cg = ftmc.CoarseGrainRNA(cg_filename)
    sd = cg.get_node_from_residue_num(stem_def.define[0])
    chain = ftup.get_biggest_chain(pdb_filename)
    chain = ftup.extract_subchain_from_res_list(chain, 
                                       list(cg.define_residue_num_iterator(sd)))

    align_chain_to_stem(cg, chain, sd, stem, use_average_method)

    for i in range(stem_def.bp_length):
        #print "i:", i
        if cg_orig.seq_ids[orig_def[0] + i - 1] in new_chain:
            new_chain.detach_child(new_chain[orig_def[0] + i].id)

        e = chain[cg.seq_ids[stem_def.define[0] + i-1]]
        e.id = cg_orig.seq_ids[orig_def[0] + i - 1]
        print "adding:", e.id
        new_chain.add(e)

        if cg_orig.seq_ids[orig_def[2] + i - 1] in new_chain:
            new_chain.detach_child(new_chain[orig_def[2] + i].id)

        e = chain[cg.seq_ids[stem_def.define[2] + i - 1]]
        e.id = cg_orig.seq_ids[orig_def[2] + i-1] #(e.id[0], orig_def[2] + i, e.id[2])
        print "adding:", e.id
        new_chain.add(e)

    return new_chain

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


def reconstruct_stem(sm, stem_name, new_chain, stem_library=dict(), stem=None):
    if stem is None:
        stem = sm.stems[stem_name]

    stem_def = sm.stem_defs[stem_name]
    orig_def = sm.bg.defines[stem_name]

    return reconstruct_stem_core(sm.bg, stem_def, orig_def, new_chain, stem_library, stem)

def place_new_stem(prev_stem, stem_params, bulge_params, (s1b, s1e), stem_name=''):
    '''
    Place a new stem with a particular orientation with respect
    to the previous stem.

    @param prev_stem: The already existing stem
    @param stem_params: A description of the new stem (StemStat)
    @param bulge_params: The AngleStat that specifies how the new stem will
                         be oriented with respect to the first stem.
    @param (s1b, s1e): Which side of the first stem to place the second on.
    '''
    stem = StemModel()
    
    stem1_basis = cuv.create_orthonormal_basis(prev_stem.vec((s1b, s1e)), prev_stem.twists[s1e]).transpose()
    start_location = cgg.stem2_pos_from_stem1_1(stem1_basis, bulge_params.position_params())
    stem_orientation = cgg.stem2_orient_from_stem1_1(stem1_basis, [stem_params.phys_length] + list(bulge_params.orientation_params()))
    twist1 = cgg.twist2_orient_from_stem1_1(stem1_basis, bulge_params.twist_params())

    mid1 = prev_stem.mids[s1e] + start_location
    mid2 = mid1 + stem_orientation

    stem.mids = (mid1, mid2)

    twist2 = cgg.twist2_from_twist1(stem_orientation, twist1, stem_params.twist_angle)
    stem.twists = (twist1, twist2)

    return stem

class SpatialModel:
    '''
    A way of building RNA structures given angle statistics as well
    as length statistics.
    '''

    def __init__(self, bg, stats_file=cbc.Configuration.stats_file, angle_defs = None, stem_defs = None, loop_defs = None):
        '''
        Initialize the structure.

        @param bg: The BulgeGraph containing information about the coarse grained
                   elements and how they are connected.
        @param angle_stats: The statistics about the inter-helical angles.
        @param angle_defs: Pre-determined statistics for each bulge
        '''

        self.stems = dict()
        self.bulges = dict()
        self.chain = bpdb.Chain.Chain(' ')
        self.build_chain = False
        self.constraint_energy = None
        self.junction_constraint_energy = None

        self.bg = bg
        self.add_to_skip()
        
        for s in bg.stem_iterator():
            try:
                cgg.add_virtual_residues(self.bg,s)
            except KeyError:
                # The structure is probably new and doesnt have coordinates yet
                continue

    def sample_stats(self):
        self.sample_angles()
        self.sample_stems()
        self.sample_loops()
        self.sample_fiveprime()
        self.sample_threeprime()

    def resample(self, d):
        if d[0] == 's':
            self.sample_stem(d)
        else:
            if len(self.bg.edges[d]) == 2:
                self.sample_angle(d)

    def sample_angle(self, d):
        size = self.bg.get_bulge_dimensions(d)

        #HACK to overcome some sparse statistics
        #TODO: remove this
        '''
        if size == (1,7):
            size = (2,7)
        if size == (2,7):
            size = (3,7)
        '''

        angle_defs = self.angle_defs

        connections = list(self.bg.edges[d])
        (s1b, s1e) = self.bg.get_sides(connections[0], d)
        (s2b, s2e) = self.bg.get_sides(connections[1], d)
        dir1 = self.bg.get_stem_direction(connections[0], connections[1])
        dir2 = self.bg.get_stem_direction(connections[1], connections[0])

        ang_type1 = self.bg.connection_type(d, connections)
        ang_type3 = -ang_type1

        try:
            angle_defs[d][ang_type1] = choice(cbs.get_angle_stats()[(size[0], size[1], ang_type1)])
            angle_defs[d][ang_type3] = choice(cbs.get_angle_stats()[(size[0], size[1], ang_type3)])
        except IndexError:
            #print >>sys.stderr, "No statistics for bulge %s of size: %s" % (d, size)

            (dist, size1, size2, _) = cbs.get_angle_stat_dims(size[0], size[1], ang_type1)[0]
            angle_defs[d][ang_type1] = choice(cbs.get_angle_stats()[(size1, size2, ang_type1)])
            #print >>sys.stderr, "Using size instead:", size1, size2
            (dist, size1, size2, _) = cbs.get_angle_stat_dims(size[0], size[1], ang_type1)[0]
            angle_defs[d][ang_type3] = choice(cbs.get_angle_stats()[(size1, size2, ang_type3)])
            #print >>sys.stderr, "Using size instead:", size1, size2
            '''
            print cbs.get_angle_stat_dims(size[0], size[1], ang_type1)
            print cbs.get_angle_stat_dims(size[0], size[1], ang_type3)
            '''

    def sample_angles(self):
        '''
        Sample statistics for each bulge region. In this case they'll be random.

        @return: A dictionary where the key is the name of a bulge and the value is a 
            statistic for the angles of that bulge.
        '''
        self.angle_defs = c.defaultdict(lambda: c.defaultdict(dict))
        '''
        angle_defs['start'][0][1] = cbs.AngleStat()
        angle_defs['start'][1][0] = cbs.AngleStat()
        '''

        for d in self.bg.defines.keys():
            if d[0] != 's': 
                if len(self.bg.edges[d]) == 2:
                    self.sample_angle(d)
                '''
                else:
                    angle_defs[d][0][0] = cbs.AngleStat()
                    pass
                '''

        #self.angle_defs = angle_defs

    def sample_stem(self, d):
        stem_defs = self.stem_defs
        define = self.bg.defines[d]
        length = self.bg.stem_length(d)

        # retrieve a random entry from the StemStatsDict collection
        ss = choice(cbs.get_stem_stats()[length])
        stem_defs[d] = ss

    def sample_stems(self):
        '''
        Sample statistics for each stem region.

        @return: A dictionary containing statistics about the stems in the structure.
        '''
        self.stem_defs = dict()

        for d in self.bg.defines.keys():
            if d[0] == 's':
                self.sample_stem(d)

    def sample_native_stems(self):
        '''
        Sample the native stems for each stem region.

        @return: A dictionary containing the real statistics about each stem in the
        structure.
        '''
        stem_defs = dict()

        for stats in cbs.get_stem_stats().values():
            for stat in stats:
                for d in self.bg.sampled.keys():
                    if d[0] == 's':
                        if stat.pdb_name == self.bg.sampled[d][0]:
                            #define = " ".join(map(str,self.bg.defines[d]))
                            if stat.define == self.bg.sampled[d][1:]:
                                stem_defs[d] = stat

        self.stem_defs = stem_defs

    def sampled_from_bg(self):
        '''
        Get the information about the sampled elements from the underlying BulgeGraph.
        '''
        # get the stem defs
        self.get_sampled_bulges()

        self.stem_defs = dict()
        for s in self.bg.stem_iterator():
            sl = self.bg.stem_length(s)
            for ss in cbs.get_stem_stats()[sl]:
                if ss.pdb_name == self.bg.sampled[s][0] and ss.define == self.bg.sampled[s][1:]:
                    self.stem_defs[s] = ss

        self.angle_defs = c.defaultdict(lambda: c.defaultdict(dict))
        for b in it.chain(self.bg.iloop_iterator(), self.bg.mloop_iterator()):
            if b not in self.bg.sampled.keys():
                print >>sys.stderr, "%s not sampled" % (b)
                continue

            size = self.bg.get_bulge_dimensions(b)
            sb = self.bg.sampled[b]

            (dist, size1, size2, _) = cbs.get_angle_stat_dims(size[0], size[1], sb[1])[0]
            size = (size1, size2)

            for ang_s in cbs.get_angle_stats()[(size[0], size[1], sb[1])]:
                if ang_s.pdb_name == sb[0] and ang_s.define == sb[2:]:
                    print >>sys.stderr, "some stuff", b

                    self.angle_defs[b][sb[1]] = ang_s

        self.loop_defs = dict()
        for l in self.bg.hloop_iterator():
            sl = self.bg.get_length(l)
            for ls in cbs.get_loop_stats()[sl]:
                if ls.pdb_name == self.bg.sampled[l][0] and ls.define == self.bg.sampled[l][1:]:
                    self.loop_defs[l] = ls

        self.fiveprime_defs = dict()
        for l in self.bg.floop_iterator():
            sl = self.bg.get_length(l)
            for ls in cbs.get_fiveprime_stats()[sl]:
                if ls.pdb_name == self.bg.sampled[l][0] and ls.define == self.bg.sampled[l][1:]:
                    self.fiveprime_defs[l] = ls

        self.threeprime_defs = dict()
        for l in self.bg.tloop_iterator():
            sl = self.bg.get_length(l)
            for ls in cbs.get_threeprime_stats()[sl]:
                if ls.pdb_name == self.bg.sampled[l][0] and ls.define == self.bg.sampled[l][1:]:
                    self.threeprime_defs[l] = ls

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

    def sample_loops(self):
        '''
        Sample statistics for each loop region.

        @return: A dictionary containing statistics about the loops in the structure.
        '''
        loop_defs = dict()

        for d in self.bg.hloop_iterator():
            define = self.bg.defines[d]
            length = self.bg.get_length(d)

            # retrieve a random entry from the StemStatsDict collection
            try:
                ls = choice(cbs.get_loop_stats()[length])
            except IndexError:
                print >>sys.stderr, "Error sampling loop %s of size %s. No available statistics." % (d, str(length))
                sys.exit(1)

            loop_defs[d] = ls

        self.loop_defs = loop_defs

    def sample_fiveprime(self):
        '''
        Sample statistics for the 5' unpaired region.

        @return: A dictionary containing statistics about the 5' unpaired regions in the structure.
        '''
        fiveprime_defs = dict()

        for d in self.bg.floop_iterator():
            define = self.bg.defines[d]
            length = self.bg.get_length(d)

            # retrieve a random entry from the StemStatsDict collection
            try:
                ls = choice(cbs.get_fiveprime_stats()[length])
            except IndexError:
                print >>sys.stderr, "Error sampling 5' %s of size %s. No available statistics." % (d, str(length))
                continue

            fiveprime_defs[d] = ls

        self.fiveprime_defs = fiveprime_defs

    def sample_threeprime(self):
        '''
        Sample statistics for the 3' unpaired region.

        @return: A dictionary containing statistics about the 3' unpaired regions in the structure.
        '''
        threeprime_defs = dict()

        for d in self.bg.tloop_iterator():
            define = self.bg.defines[d]
            length = self.bg.get_length(d)

            # retrieve a random entry from the StemStatsDict collection
            try:
                ls = choice(cbs.get_threeprime_stats()[length])
            except IndexError:
                print >>sys.stderr, "Error sampling threeprime %s of size %s. No available statistics." % (d, str(length))
                continue

            threeprime_defs[d] = ls

        self.threeprime_defs = threeprime_defs

    def add_loop(self, name, prev_stem_node, params=None, loop_defs=None):
        '''
        Connect a loop to the previous stem.
        '''
        if loop_defs == None:
            loop_defs = self.loop_defs

        prev_stem = self.stems[prev_stem_node]
        (s1b, s1e) = self.bg.get_sides(prev_stem_node, name)

        if params == None:
            r = loop_defs[name].phys_length
            u = loop_defs[name].u
            v = loop_defs[name].v

            params = (r, u, v)

        start_mid = prev_stem.mids[s1b]
        (r, u, v) = params

        direction = cgg.stem2_pos_from_stem1(prev_stem.vec((s1e, s1b)), prev_stem.twists[s1b], (r, u, v))
        end_mid = start_mid + direction
        self.bulges[name] = BulgeModel((start_mid, end_mid))

    def find_start_node(self):
        '''
        Find a node from which to begin building. This should ideally be a loop
        region as from there we can just randomly orient the first stem.
        '''

        edge = self.bg.sorted_stem_iterator().next()
        define = 'start'
        return (edge, define, StemModel(edge))


        for define in self.bg.defines.keys():
            if define[0] == 'h' or define[0] == 'f' or define[0] == 't':
                for edge in self.bg.edges[define]:
                    return (edge, define, StemModel(edge))

    def save_sampled_stems(self):
        '''
        Save the information about the sampled stems to the bulge graph file.
        '''
        for sd in self.stem_defs.items():
            self.bg.sampled[sd[0]] = [sd[1].pdb_name] + sd[1].define

    def save_sampled_angles(self):
        for (bulge, ang_type) in self.sampled_bulge_sides:
            ad = self.angle_defs[bulge][ang_type]
            self.bg.sampled[bulge] = [ad.pdb_name] + [ang_type] + ad.define

    def save_sampled_loops(self):
        for sd in self.loop_defs.items():
            self.bg.sampled[sd[0]] = [sd[1].pdb_name] + sd[1].define

    def save_sampled_fiveprimes(self):
        for sd in self.fiveprime_defs.items():
            self.bg.sampled[sd[0]] = [sd[1].pdb_name] + sd[1].define

    def save_sampled_threeprimes(self):
        for sd in self.threeprime_defs.items():
            self.bg.sampled[sd[0]] = [sd[1].pdb_name] + sd[1].define

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

        mat = cuv.get_double_alignment_matrix(target, [vec1, twist1])

        return mat

    def get_random_stem_stats(self, name):
        '''
        Return a random set of parameters with which to create a stem.
        '''

        return self.stem_defs[name]

    def get_random_bulge_stats(self, name, ang_type):
        '''
        Return a random set of parameters with which to create a bulge.
        '''
        #if name[0] != 's' and self.bg.weights[name] == 1 and len(self.bg.edges[name]) == 1:
        if name[0] == 'h':
            return cbs.AngleStat()

        return self.angle_defs[name][ang_type]

    def add_stem(self, stem_name, stem_params, prev_stem, bulge_params, (s1b, s1e)):
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

        stem = place_new_stem(prev_stem, stem_params, bulge_params, (s1b, s1e), stem_name)

        stem.name = stem_name

        if self.build_chain:
            reconstruct_stem(self, stem_name, self.chain, stem_library=cbc.Configuration.stem_library, stem=stem)

        return stem

    def fill_in_bulges_and_loops(self):
        loops = list(self.bg.hloop_iterator())
        fiveprime = list(self.bg.floop_iterator())
        threeprime = list(self.bg.tloop_iterator())

        for d in self.bg.defines.keys():
            if d[0] != 's':
                if d in loops:
                    self.add_loop(d, list(self.bg.edges[d])[0])
                     #add loop
                    pass
                elif d in fiveprime:
                    self.add_loop(d, list(self.bg.edges[d])[0],
                                  loop_defs = self.fiveprime_defs)
                elif d in threeprime:
                    self.add_loop(d, list(self.bg.edges[d])[0],
                                  loop_defs = self.threeprime_defs)
                else:
                    connections = list(self.bg.edges[d])

                    # Should be a bulge connecting two stems
                    assert(len(connections) == 2)
                    for conn in connections:
                        assert(conn[0] == 's')

                    (s1b, s1e) = self.bg.get_sides(connections[0], d)
                    (s2b, s2e) = self.bg.get_sides(connections[1], d)

                    s1mid = self.stems[connections[0]].mids[s1b]
                    s2mid = self.stems[connections[1]].mids[s2b]

                    self.bulges[d] = BulgeModel((s1mid, s2mid))
                    self.closed_bulges += [d]

    def stem_to_coords(self, stem):
        sm = self.stems[stem]

        self.bg.coords[stem] = (sm.mids[0], sm.mids[1])
        self.bg.twists[stem] = (sm.twists[0], sm.twists[1])

        '''
        for edge in self.bg.edges[stem]:
            if self.bg.weights[edge] == 2:
                cgg.add_virtual_residues(self.bg, edge)
        '''

        cgg.add_virtual_residues(self.bg, stem)

    def elements_to_coords(self):
        '''
        Add all of the stem and bulge coordinates to the BulgeGraph data structure.
        '''

        #for stem in self.stems.keys():
        for stem in self.newly_added_stems:
            self.stem_to_coords(stem)

        for bulge in self.bulges.keys():
            bm = self.bulges[bulge]

            self.bg.coords[bulge] = (bm.mids[0], bm.mids[1])

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

    def finish_building(self):
        self.fill_in_bulges_and_loops()
        self.elements_to_coords()
        self.save_sampled_stems()
        self.save_sampled_angles()
        self.save_sampled_loops()
        self.save_sampled_fiveprimes()
        self.save_sampled_threeprimes()

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

    def traverse_and_build(self, start='start'):
        '''
        Build a 3D structure from the graph in self.bg.

        This is done by doing a breadth-first search through the graph
        and adding stems. 

        Once all of the stems have been added, the bulges and loops
        are added.

        If a constraint energy is provided, then the nascent structure
        must fullfill the constraint at every step of the process.
        '''
        constraint_energy = self.constraint_energy
        '''
        import traceback
        print "".join(traceback.format_stack()[-3:])
        '''

        #print >>sys.stderr, "traverse_and_build"
        self.visited = set()
        self.to_visit = []
        #self.stems = dict()
        #self.bulges = dict()
        self.sampled_bulges = []
        self.sampled_bulge_sides = []
        self.closed_bulges = []
        self.newly_added_stems = []

        new_visited = []
        # the start node should be a loop region
        self.to_visit = [self.find_start_node()]
        paths = c.defaultdict(list)

        self.visit_order = []
        self.prev_stem_list = []

        counter = 0
        '''
        self.bg.coords = dict()
        self.bg.bases = dict()
        self.bg.stem_invs = dict()
        '''
        started = False

        if start == '':
            started = True

        restart=False

        #print "starting:"

        while True:
            while len(self.to_visit) > 0:
                #self.to_visit.sort(key=lambda x: ('a' if x[1][0] == 's' else x[1][0], -self.bg.stem_length(x[1])))
                (curr_node, prev_node, prev_stem) = self.to_visit.pop(0)
                tbc = True
                while curr_node in self.visited or curr_node in self.to_skip:
                    if len(self.to_visit) > 0:
                        (curr_node, prev_node, prev_stem) = self.to_visit.pop()
                    else:
                        #self.finish_building()
                        tbc = False
                        break

                if not tbc:
                    break

                # keep track of the leaf to root paths
                paths[curr_node] += [curr_node]
                paths[curr_node] += paths[prev_node]

                self.visited.add(curr_node)

                v = list(self.visited)
                v.sort()

                stem = prev_stem

                #if curr_node == start:
                if start in paths[curr_node] or start == 'start':
                    started = True

                if curr_node[0] == 's':
                    params = self.get_random_stem_stats(curr_node)

                    if prev_node == 'start':
                        (s1b, s1e) = (0, 1)
                    else:
                        (s1b, s1e) = self.bg.get_sides(curr_node, prev_node)

                    # get some parameters for the previous bulge
                    if prev_node == 'start':
                        (ps1b, ps1e) = (1, 0)
                        ang_type = 1
                        prev_params = cbs.AngleStat()
                    else:
                        (ps1b, ps1e) = self.bg.get_sides(prev_stem.name, prev_node)
                        ang_type = self.bg.connection_type(prev_node, 
                                                           [prev_stem.name, 
                                                            curr_node])

                        prev_params = self.get_random_bulge_stats(prev_node,
                                                                 ang_type)

                    self.sampled_bulges += [prev_node]
                    if len(self.bg.edges[prev_node]) == 2:
                        self.sampled_bulge_sides += [(prev_node, ang_type)]

                    # the previous stem should always be in the direction(0, 1) 
                    if started:
                        new_visited += [curr_node]
                        '''
                        if curr_node == 's10':
                            print "started prev_node:", prev_node, "curr_node:", curr_node, "start:", start
                        '''
                        #stem = self.add_stem(curr_node, params, prev_stem, prev_params, (0, 1))
                        #print "ps1b:", ps1b, "ps1e", ps1e
                        self.visit_order += [prev_node]

                        self.prev_stem_list += [prev_stem.name]
                        stem = self.add_stem(curr_node, params, prev_stem, prev_params, (ps1e, ps1b))
                        self.newly_added_stems += [curr_node]

                        # the following is done to maintain the invariant that mids[s1b] is
                        # always in the direction of the bulge from which s1b was obtained
                        # i.e. s1e -> s1b -> bulge -> s2b -> s2e
                        #self.stems[curr_node] = stem

                        if s1b == 1:
                            self.stems[curr_node] = stem.reverse()
                        else:
                            self.stems[curr_node] = stem


                        if constraint_energy != None and not restart:
                            self.stem_to_coords(curr_node)
                            e1 = constraint_energy.eval_energy(self, nodes=self.visited, new_nodes = new_visited)
                            #e1 = constraint_energy.eval_energy(self)
                            if e1 > 10:
                                #self.bg.to_file('bad1.cg')
                                #print >>sys.stderr, "exiting0", e1
                                #sys.exit(1)

                                bb = set(self.constraint_energy.bad_bulges)
                                bp = []
                                for b in bb:
                                    bp += [(len(paths[b]), b)]

                                bp.sort()
                                sb = max((bp))[1]
                                #sb = random.choice(list(bb))

                                for p in paths[sb]:
                                    if p[0] == 's':
                                        continue
                                    if random.random() < 0.8:
                                        break

                                to_change = p
                                restart = True
                                break
                            else:
                                pass
                            new_visited = []

                    else:
                        '''
                        if curr_node == 's13':
                            print "unstarted prev_node:", prev_node, "start:", start
                        '''
                        if s1b == 1:
                            stem = self.stems[curr_node].reverse()
                        else:
                            stem = self.stems[curr_node]

                to_append = []
                for edge in self.bg.edges[curr_node]:
                    if edge not in self.visited:
                        to_append.append((edge, curr_node, stem))
                        to_append.sort(key=lambda x: -self.bg.stem_length(x[1]))

                self.to_visit += to_append

                counter += 1

            if not restart and self.constraint_energy != None:
                e1 = self.constraint_energy.eval_energy(self, nodes=self.visited, new_nodes = None)
                if e1 > 0.:
                    #self.bg.to_file('bad.cg')
                    #print >>sys.stderr, "exiting1", e1
                    #sys.exit(1)
                    bb = set(self.constraint_energy.bad_bulges)
                    bp = []
                    for b in bb:
                        bp += [(len(paths[b]), b)]

                    bp.sort()
                    sb = max((bp))[1]
                    #sb = random.choice(list(bb))

                    for p in paths[sb]:
                        if p[0] == 's':
                            continue
                        if random.random() < 0.5:
                            break

                    to_change = p
                    restart = True

            if not restart and self.junction_constraint_energy != None:
                e1 = self.junction_constraint_energy.eval_energy(self)
                if e1 > 0.:
                    #self.bg.to_file('bad2.cg')
                    #print >>sys.stderr, "exiting2", e1
                    #sys.exit(1)
                    to_change = random.choice(self.junction_constraint_energy.bad_bulges)
                    restart = True

            if restart:
                self.resample(to_change)
                self.sampled_bulges = []
                self.sampled_bulge_sides = []
                self.closed_bulges = []
                self.newly_added_stems = []
                self.visited = set()
                self.to_visit = [self.find_start_node()]
                self.visit_order = []
                #paths = c.defaultdict(list)
                paths = c.defaultdict(list)
                start = to_change
                #start = 'start'
                started = False
                new_visited = []
                #sys.exit(1)
                #self.traverse_and_build(to_change)
                #return
                restart = False
            else:
                break

        self.finish_building()

