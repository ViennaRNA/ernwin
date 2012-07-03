#!/usr/bin/python

from corgy.utilities.data_structures import DefaultDict
from corgy.visual.pymol import PymolPrinter
from corgy.builder.stats import AngleStat, AngleStatsDict
from corgy.builder.stats import StemStat, StemStatsDict
from corgy.builder.stats import get_loop_length
from corgy.graph.graph_pdb import stem2_pos_from_stem1
from corgy.graph.graph_pdb import stem2_orient_from_stem1
from corgy.graph.graph_pdb import twist2_orient_from_stem1
from corgy.graph.graph_pdb import twist2_from_twist1
from corgy.utilities.vector import get_double_alignment_matrix

from numpy import array
from random import choice, uniform
from sys import stderr
from math import pi


class StemModel:
    '''
    A way of encapsulating the coarse grain 3D stem.
    '''

    def __init__(self, mids=None, twists=None):
        if mids == None:
            self.mids = (array([0., 0., 0.]), array([0., 0., 1.0]))
        else:
            self.mids = mids
        if twists == None:
            self.twists = (array([0., 1., 0.]), array([1., 0., 0.0]))
        else:
            self.twists = twists

    def __str__(self):
        return str(self.mids) + '\n' + str(self.twists)

    def reverse(self):
        '''
        Reverse this stem's orientation so that the order of the mids
        is backwards. I.e. mids[1] = mids[0]...
        '''
        return StemModel((self.mids[1], self.mids[0]), (self.twists[1], self.twists[0]))

    def vec(self, (from_side, to_side) = (0, 1)):
        return self.mids[to_side] - self.mids[from_side]

class BulgeModel:
    '''
    A way of encapsulating a coarse grain 3D loop.
    '''

    def __init__(self, mids=None):
        if mids == None:
            self.mids = (array([0., 0., 0.]), array([0., 0., 1.0]))
        else:
            self.mids = mids

    def __str__(self):
        return str(self.mids)
            

class SpatialModel:
    '''
    A way of building RNA structures given angle statistics as well
    as length statistics.
    '''

    def __init__(self, bg, angle_stats, stem_stats):
        '''
        Initialize the structure.

        @param bg: The BulgeGraph containing information about the coarse grained
                   elements and how they are connected.
        @param angle_stats: The statistics about the inter-helical angles.
        '''
        self.angle_stats = angle_stats
        self.stem_stats = stem_stats

        self.bg = bg

        self.pymol_printer = PymolPrinter()

        pass

    def add_loop(self, name, prev_stem_node, params=None):
        '''
        Connect a loop to the previous stem.
        '''
        prev_stem = self.stems[prev_stem_node]
        (s1b, s1e) = self.bg.get_sides(prev_stem_node, name)

        if params == None:
            length = get_loop_length(self.bg, name)
            u = uniform(0, pi)
            v = uniform(-pi/2, pi/2)

            params = (length, u, v)

        start_mid = prev_stem.mids[s1b]
        (r, u, v) = params

        direction = stem2_pos_from_stem1(prev_stem.vec((s1e, s1b)), prev_stem.twists[s1b], (r, u, v))
        end_mid = start_mid + direction
        self.bulges[name] = BulgeModel((start_mid, end_mid))


    def find_start_node(self):
        '''
        Find a node from which to begin building. This should ideally be a loop
        region as from there we can just randomly orient the first stem.
        '''
        for define in self.bg.defines.keys():
            if define[0] != 's' and self.bg.weights[define] == 1 and len(self.bg.edges[define]) == 1:

                length = get_loop_length(self.bg, define)
                params = (length, 0., 0.)

                prev_stem = StemModel()

                #self.add_loop(define, StemModel(), params)

                ang_stats = AngleStat()
                ang_stats.pdb_name = 'start'
                ang_stats.r1 = length

                for edge in self.bg.edges[define]:
                    self.to_visit += [(edge, define, StemModel(), ang_stats)] 

                break


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

        target = [array([1., 0., 0.]), array([0., 1., 0.])]


        vec1 = self.stems[edge].vec()
        twist1 = self.stems[edge].twists[1]

        mat = get_double_alignment_matrix(target, [vec1, twist1])

        return mat

    def get_random_stem_stats(self, name):
        '''
        Return a random set of parameters with which to create a stem.
        '''
        # get the stem's length

        d = self.bg.defines[name]
        length = abs(d[1] - d[0])

        # retrieve a random entry from the StemStatsDict collection
        ss = choice(self.stem_stats[length])

        return ss

    def get_random_bulge_stats(self, name):
        '''
        Return a random set of parameters with which to create a bulge.
        '''

        # get the bulge's size
        size = self.bg.get_bulge_dimensions(name)

        # retrieve a random entry from the AngleStatsDict collection
        stats = choice(self.angle_stats[size[0]][size[1]])
        return stats

    def add_stem(self, stem_params, prev_stem, bulge_params, (s1b, s1e)):
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
        start_location = stem2_pos_from_stem1(prev_stem.vec((s1b, s1e)), prev_stem.twists[s1e], bulge_params.position_params())
        stem_orientation = stem2_orient_from_stem1(prev_stem.vec((s1b, s1e)), prev_stem.twists[s1e], [stem_params.phys_length] + list(bulge_params.orientation_params()))
        twist1 = twist2_orient_from_stem1(prev_stem.vec((s1b, s1e)), prev_stem.twists[s1e], bulge_params.twist_params())
    
        stem = StemModel()

        mid1 = prev_stem.mids[s1e] + start_location
        mid2 = mid1 + stem_orientation

        stem.mids = (mid1, mid2)

        twist2 = twist2_from_twist1(stem_orientation, twist1, stem_params.twist_angle)

        stem.twists = (twist1, twist2)

        return stem

    def fill_in_bulges_and_loops(self):
        for d in self.bg.defines.keys():
            if d[0] != 's':
                if len(self.bg.edges[d]) == 1:
                    self.add_loop(d, list(self.bg.edges[d])[0])
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

    def elements_to_coords(self):
        '''
        Add all of the stem and bulge coordinates to the BulgeGraph data structure.
        '''

        for stem in self.stems.keys():
            sm = self.stems[stem]

            self.bg.coords[stem] = (sm.mids[0], sm.mids[1])
            self.bg.twists[stem] = (sm.twists[0], sm.twists[1])

        for bulge in self.bulges.keys():
            bm = self.bulges[bulge]

            self.bg.coords[bulge] = (bm.mids[0], bm.mids[1])

    def traverse_and_build(self):
        '''
        Build a 3D structure from the graph in self.bg.

        This is done by doing a breadth-first search through the graph
        and adding stems. 

        Once all of the stems have been added, the bulges and loops
        are added.
        '''

        self.visited = set()
        self.to_visit = []
        self.stems = dict()
        self.bulges = dict()

        # the start node should be a loop region
        self.find_start_node()

        counter = 0
        self.bg.coords = dict()

        while len(self.to_visit) > 0:
            (curr_node, prev_node, prev_stem, prev_params) = self.to_visit.pop(0)

            while curr_node in self.visited:
                if len(self.to_visit) > 0:
                    (curr_node, prev_node, prev_stem, prev_params) = self.to_visit.pop(0)
                else:
                    return
                    
            self.visited.add(curr_node)
            stem = prev_stem

            if curr_node[0] == 's':
                params = self.get_random_stem_stats(curr_node)
                (s1b, s1e) = self.bg.get_sides(curr_node, prev_node)

                # the previous stem should always be in the direction(0, 1) 
                stem = self.add_stem(params, prev_stem, prev_params, (0, 1))

                # the following is done to maintain the invariant that mids[s1b] is
                # always in the direction of the bulge from which s1b was obtained
                # i.e. s1e -> s1b -> bulge -> s2b -> s2e
                if s1b == 1:
                    self.stems[curr_node] = stem.reverse()
                else:
                    self.stems[curr_node] = stem

            else:
                if len(self.bg.edges[curr_node]) > 1:
                    params = self.get_random_bulge_stats(curr_node)

            for edge in self.bg.edges[curr_node]:
                if edge in self.visited and edge != prev_node and edge[0] == 's':
                    next_stem = self.stems[edge]

                if edge not in self.visited:
                    self.to_visit.append((edge, curr_node, stem, params))

            counter += 1

        self.fill_in_bulges_and_loops()
        
        self.elements_to_coords()


