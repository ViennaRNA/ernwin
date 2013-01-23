#!/usr/bin/python

import Bio.PDB as bpdb
import os
import numpy as np
import numpy.linalg as nl
import math
import sys
import collections as c

import corgy.builder.config as cbc
import corgy.builder.stats as cbs
import corgy.graph.graph_pdb as cgg
import corgy.utilities.vector as cuv
import corgy.utilities.debug as cud

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

def define_to_stem_model(chain, define):
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

    mids = cgg.get_mids(chain, define)

    stem.mids = tuple([m.get_array() for m in mids])
    stem.twists = cgg.get_twists(chain, define)

    return stem

def get_stem_rotation_matrix(stem, (u, v, t)):
    twist1 = stem.twists[0]

    # rotate around the stem axis to adjust the twist

    # rotate down from the twist axis
    comp1 = np.cross(stem.vec(), twist1)

    rot_mat1 = cuv.rotation_matrix(stem.vec(), t)
    rot_mat2 = cuv.rotation_matrix(twist1, u - math.pi/2)
    rot_mat3 = cuv.rotation_matrix(comp1, v)

    rot_mat4 = np.dot(rot_mat3, np.dot(rot_mat2, rot_mat1))

    return rot_mat4

def align_chain_to_stem(chain, define, stem2):
    stem1 = define_to_stem_model(chain, define)

    (r, u, v, t) = cgg.get_stem_orientation_parameters(stem1.vec(), stem1.twists[0], stem2.vec(), stem2.twists[0])
    rot_mat = get_stem_rotation_matrix(stem1, (math.pi-u, -v, -t))
    rotate_chain(chain, np.linalg.inv(rot_mat), stem1.mids[0])
    translate_chain(chain, stem2.mids[0] - stem1.mids[0])

def reconstruct_stem(sm, stem_name, new_chain, stem_library=dict(), stem=None):
    '''
    Reconstruct a particular stem.
    '''
    if stem is None:
        stem = sm.stems[stem_name]

    stem_def = sm.stem_defs[stem_name]

    filename = '%s_%s.pdb' % (stem_def.pdb_name, "_".join(map(str, stem_def.define)))
    cud.pv('filename')
    #print "stem_name:", stem_name, "stem_def:", stem_def, "filename:", filename
    pdb_file = os.path.join(cbc.Configuration.stem_fragment_dir, filename)

    #print len(stem_library.keys())
    if filename in stem_library.keys():
    #if False:
        chain = stem_library[filename].copy()
    else:
        chain = list(bpdb.PDBParser().get_structure('temp', pdb_file).get_chains())[0]
        stem_library[filename] = chain.copy()

    align_chain_to_stem(chain, stem_def.define, stem)

    for i in range(stem_def.bp_length+1):
        #print "i:", i
        if sm.bg.defines[stem_name][0] + i in new_chain:
            new_chain.detach_child(new_chain[sm.bg.defines[stem_name][0] + i].id)

        e = chain[stem_def.define[0] + i]
        e.id = (e.id[0], sm.bg.defines[stem_name][0] + i, e.id[2])
        #print "adding:", e.id
        new_chain.add(e)

        if sm.bg.defines[stem_name][2] + i in new_chain:
            new_chain.detach_child(new_chain[sm.bg.defines[stem_name][2] + i].id)

        e = chain[stem_def.define[2] + i]
        e.id = (e.id[0], sm.bg.defines[stem_name][2] + i, e.id[2])
        #print "adding:", e.id
        new_chain.add(e)

def place_new_stem(prev_stem, stem_params, bulge_params, (s1b, s1e)):
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

        self.angle_stats = cbs.get_angle_stats()
        self.stem_stats = cbs.get_stem_stats()
        self.loop_stats = cbs.get_loop_stats()
        self.stems = dict()
        self.bulges = dict()
        self.chain = bpdb.Chain.Chain(' ')
        self.build_chain = False

        self.bg = bg

        if self.angle_stats == None:
            return

    def sample_stats(self):
        self.sample_angles()
        self.sample_stems()
        self.sample_loops()

    def sample_angles(self):
        '''
        Sample statistics for each bulge region. In this case they'll be random.

        @return: A dictionary where the key is the name of a bulge and the value is a 
            statistic for the angles of that bulge.
        '''
        angle_defs = c.defaultdict(lambda: c.defaultdict(dict))
        angle_defs['start'][0][1] = cbs.AngleStat()
        angle_defs['start'][1][0] = cbs.AngleStat()

        for d in self.bg.defines.keys():
            if d[0] != 's': 
                if len(self.bg.edges[d]) == 2:
                    size = self.bg.get_bulge_dimensions(d)

                    #HACK to overcome some sparse statistics
                    #TODO: remove this
                    '''
                    if size == (1,7):
                        size = (2,7)
                    if size == (2,7):
                        size = (3,7)
                    '''

                    connections = list(self.bg.edges[d])
                    (s1b, s1e) = self.bg.get_sides(connections[0], d)
                    (s2b, s2e) = self.bg.get_sides(connections[1], d)

                    try:
                        angle_defs[d][s1b][s2b] = choice(self.angle_stats[size[0]][size[1]][s1b][s2b])
                        angle_defs[d][s2b][s1b] = choice(self.angle_stats[size[0]][size[1]][s2b][s1b])
                    except IndexError:
                        print >>sys.stderr, "No statistics for bulge %s of size: %s" % (d, size)
                else:
                    angle_defs[d][0][0] = cbs.AngleStat()
                    pass

        self.angle_defs = angle_defs

    def sample_stems(self):
        '''
        Sample statistics for each stem region.

        @return: A dictionary containing statistics about the stems in the structure.
        '''
        stem_defs = dict()

        for d in self.bg.defines.keys():
            if d[0] == 's':
                define = self.bg.defines[d]
                length = abs(define[1] - define[0])

                # retrieve a random entry from the StemStatsDict collection
                ss = choice(self.stem_stats[length])
                stem_defs[d] = ss

        self.stem_defs = stem_defs

    def sample_native_stems(self):
        '''
        Sample the native stems for each stem region.

        @return: A dictionary containing the real statistics about each stem in the
        structure.
        '''
        stem_defs = dict()

        '''
        for d in self.bg.sampled_stems.keys():
            stem_defs[d] = self.bg.get_stem_stats(d) 

        '''
        for stats in self.stem_stats.values():
            for stat in stats:
                for d in self.bg.sampled_stems.keys():
                    if d[0] == 's':
                        if stat.pdb_name == self.bg.sampled_stems[d][0]:
                            #define = " ".join(map(str,self.bg.defines[d]))
                            if stat.define == self.bg.sampled_stems[d][1:]:
                                stem_defs[d] = stat

        self.stem_defs = stem_defs

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
        return
        loop_defs = dict()

        for d in self.bg.defines.keys():
            if d[0] != 's' and len(self.bg.edges[d]) == 1:
                define = self.bg.defines[d]
                length = abs(define[1] - define[0])

                # retrieve a random entry from the StemStatsDict collection
                try:
                    ls = choice(self.loop_stats[length])
                except IndexError:
                    print >>sys.stderr, "Error sampling loop %s of size %s. No available statistics." % (d, str(length))
                    sys.exit(1)

                loop_defs[d] = ls

        self.loop_defs = loop_defs

    def add_loop(self, name, prev_stem_node, params=None):
        '''
        Connect a loop to the previous stem.
        '''
        prev_stem = self.stems[prev_stem_node]
        (s1b, s1e) = self.bg.get_sides(prev_stem_node, name)

        if params == None:
            length = self.loop_defs[name].phys_length
            u = uniform(0, pi)
            v = uniform(-pi/2, pi/2)

            params = (length, u, v)

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
        for define in self.bg.defines.keys():
            if define[0] != 's' and self.bg.weights[define] == 1 and len(self.bg.edges[define]) == 1:
                for edge in self.bg.edges[define]:
                    if self.bg.get_sides(edge, define)[0] == 0:
                        return (edge, define, StemModel(edge))


    def save_sampled_stems(self):
        '''
        Save the information about the sampled stems to the bulge graph file.
        '''
        for sd in self.stem_defs.items():
            self.bg.sampled_stems[sd[0]] = [sd[1].pdb_name] + sd[1].define
            #print "self.bg.sampled_stems[sd[0]]:", self.bg.sampled_stems[sd[0]]

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

    def get_random_bulge_stats(self, name, (s1b, s2b)):
        '''
        Return a random set of parameters with which to create a bulge.
        '''
        return self.angle_defs[name][s1b][s2b]

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
        stem = place_new_stem(prev_stem, stem_params, bulge_params, (s1b, s1e))
        stem.name = stem_name

        if self.build_chain:
            reconstruct_stem(self, stem_name, self.chain, stem_library=cbc.Configuration.stem_library, stem=stem)

        return stem

    def fill_in_bulges_and_loops(self):
        for d in self.bg.defines.keys():
            if d[0] != 's':
                if len(self.bg.edges[d]) == 1:
                    #self.add_loop(d, list(self.bg.edges[d])[0])
                    # add loop
                    pass
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

    def elements_to_coords(self):
        '''
        Add all of the stem and bulge coordinates to the BulgeGraph data structure.
        '''

        #for stem in self.stems.keys():
        for stem in self.newly_added_stems:
            sm = self.stems[stem]

            self.bg.coords[stem] = (sm.mids[0], sm.mids[1])
            self.bg.twists[stem] = (sm.twists[0], sm.twists[1])

            '''
            for edge in self.bg.edges[stem]:
                if self.bg.weights[edge] == 2:
                    cgg.add_virtual_residues(self.bg, edge)
            '''

            cgg.add_virtual_residues(self.bg, stem)

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

    def traverse_and_build(self, start=''):
        '''
        Build a 3D structure from the graph in self.bg.

        This is done by doing a breadth-first search through the graph
        and adding stems. 

        Once all of the stems have been added, the bulges and loops
        are added.
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

        # the start node should be a loop region
        self.to_visit = [self.find_start_node()]

        counter = 0
        '''
        self.bg.coords = dict()
        self.bg.bases = dict()
        self.bg.stem_invs = dict()
        '''
        started = False

        if start == '':
            started = True


        while len(self.to_visit) > 0:
            self.to_visit.sort(key=lambda x: -self.bg.stem_length(x[0]))
            (curr_node, prev_node, prev_stem) = self.to_visit.pop()

            while curr_node in self.visited:
                if len(self.to_visit) > 0:
                    (curr_node, prev_node, prev_stem) = self.to_visit.pop()
                else:
                    self.finish_building()
                    return
                    
            self.visited.add(curr_node)
            stem = prev_stem

            if prev_node == start:
                started = True

            if curr_node[0] == 's':
                params = self.get_random_stem_stats(curr_node)

                if prev_node == 'start':
                    (s1b, s1e) = (1, 0)
                else:
                    (s1b, s1e) = self.bg.get_sides(curr_node, prev_node)


                # get some parameters for the previous bulge
                #print "curr_node:", curr_node, "prev_stem.name", prev_stem.name
                (ps1b, ps1e) = self.bg.get_sides(prev_stem.name, prev_node)

                prev_params = self.get_random_bulge_stats(prev_node, (ps1b, s1b))
                self.sampled_bulges += [prev_node]
                if len(self.bg.edges[prev_node]) == 2:
                    self.sampled_bulge_sides += [(prev_node, (ps1b, s1b))]

                # the previous stem should always be in the direction(0, 1) 
                if started:
                    #stem = self.add_stem(curr_node, params, prev_stem, prev_params, (0, 1))
                    #print "ps1b:", ps1b, "ps1e", ps1e
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
                else:
                    if s1b == 1:
                        stem = self.stems[curr_node].reverse()
                    else:
                        stem = self.stems[curr_node]

            for edge in self.bg.edges[curr_node]:
                if edge not in self.visited:
                    self.to_visit.append((edge, curr_node, stem))

            counter += 1
        self.finish_building()

