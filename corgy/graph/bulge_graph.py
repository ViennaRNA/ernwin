#!/usr/bin/python

import sys, collections
import collections as c
import math
import random
import numpy as np

import corgy.builder.stats as cbs
import corgy.graph.graph_pdb as cgg
import corgy.utilities.data_structures as cuds
import corgy.utilities.vector as cuv
import corgy.utilities.debug as cud

def error_exit(message):
    print >> sys.stderr, message
    sys.exit(1)

# A wrapper for a simple dictionary addition
# Added so that debugging can be made easier
def add_bulge(bulges, bulge, context, message):
    #print >>sys.stderr,"Adding bulge", context, bulge, message
    #bulge = (context, bulge)
    bulges[context] = bulges.get(context, []) + [bulge]
    return bulges

def any_difference_of_one(stem, bulge):
    '''
    See if there's any difference of one between the two
    ends of the stem [(a,b),(c,d)] and a bulge (e,f)

    :param stem: A couple of couples (2 x 2-tuple) indicating the start and end
                 nucleotides of the stem in the form ((s1, e1), (s2, e2))
    :param bulge: A couple (2-tuple) indicating the first and last position
                  of the bulge.

    :return: True if there is an overlap between the stem nucleotides and the 
                  bulge nucleotides
             False otherwise
    '''
    for stem_part in stem:
        for part in stem_part:
            for bulge_part in bulge:
                if abs(bulge_part - part) == 0:
                    return True
    return False


def print_bulge_graph(graph):
    '''
    Print out the connections in the graph.

    :param graph: A dictionary indexed by stem number containing a set
                  of the bulge numbers that it is connected to.
    '''
    for key in graph.keys():
        stem_str = "connect s%d" % (key)
        for item in graph[key]:
            stem_str += " b%d" % (item)
        print stem_str

def print_stems(stems):
    '''
    Print the names and definitions of the stems.

    :param stems: A list of tuples of tuples of the form [((s1, e1), (s2, e2))]
                  where s1 and e1 are the nucleotides at one end of the stem
                  and s2 and e2 are the nucleotides at the other.
    '''
    for i in range(len(stems)):
        # one is added to each coordinate to make up for the fact that residues are 1-based
        ss1 = stems[i][0][0]+1
        ss2 = stems[i][0][1]+1
        se1 = stems[i][1][0]+1
        se2 = stems[i][1][1]+1

        print "define s%d 0 %d %d %d %d" % (i, min(ss1,se1), max(ss1,se1), min(ss2,se2), max(ss2,se2))

def print_bulges(bulges):
    '''
    Print the names and definitions of the bulges.

    :param bulges: A list of tuples of the form [(s, e)] where s and e are the 
                   numbers of the nucleotides at the start and end of the bulge.
    '''
    for i in range(len(bulges)):
            #print "bulge:", bulge
        bulge_str = "define b%d 1" % (i)
        bulge = bulges[i]
        #print >>sys.stderr, "bulge:", bulge
        bulge_str += " %d %d" % (bulge[0]+1, bulge[1]+1)
        print bulge_str

def condense_stem_pairs(stem_pairs):
    '''
    Given a list of stem pairs, condense them into stem definitions

    I.e. the pairs (0,10),(1,9),(2,8),(3,7) can be condensed into
    just the ends of the stem: [(0,10),(3,7)]

    :param stem_pairs: A list of tuples containing paired base numbers.

    :returns: A list of tuples of tuples of the form [((s1, e1), (s2, e2))]
                  where s1 and e1 are the nucleotides at one end of the stem
                  and s2 and e2 are the nucleotides at the other.
    '''
    stem_pairs.sort()

    prev_pair = (-10, -10)

    stems = []
    start_pair = None

    for pair in stem_pairs:
        # There's a potential bug here since we don't check the direction
        # but hopefully it won't bite us in the ass later
        if abs(pair[0] - prev_pair[0]) != 1 or abs(pair[1] - prev_pair[1]) != 1:
            if start_pair != None:
                stems += [(start_pair, prev_pair)]
            start_pair = pair
    
        prev_pair = pair

    if start_pair != None:
        stems += [(start_pair, prev_pair)]

    return stems

def print_brackets(brackets):
    '''
    Print the brackets and a numbering, for debugging purposes

    :param brackets: A string with the dotplot passed as input to this script.
    '''
    numbers = [chr(ord('0') + i % 10) for i in range(len(brackets))]
    tens = [chr(ord('0') + i / 10) for i in range(len(brackets))]
    print "brackets:\n", brackets, "\n", "".join(tens), "\n" ,"".join(numbers)

def find_bulges_and_stems(brackets):
    '''
    Iterate through the structure and enumerate the bulges and the stems that are
    present.

    The returned stems are of the form [[(s1, s2), (e1,e2)], [(s1,s2),(e1,e2)],...]
    where (s1,s2) are the residue numbers of one end of the stem and (e1,e2) are the
    residue numbers at the other end of the stem
    (see condense_stem_pairs)

    The returned bulges are of the form [(s,e), (s,e),...] where s is the start of a bulge
    and e is the end of a bulge

    :param brackets: A string with the dotplot passed as input to this script.
    '''
    prev = '.'
    context = 0

    bulges = dict()
    finished_bulges = []
    context_depths = dict()

    opens = []
    stem_pairs = []

    stems = dict()

    dots_start = 0
    dots_end = 0

    context_depths[0] = 0

    i = 0
    for i in range(len(brackets)):
        #print >> sys.stderr, "bracket:", brackets[i]
        if brackets[i] == '(':
            opens.append(i)

            if prev == '(':
                context_depths[context] = context_depths.get(context, 0) + 1
                continue
            else:
                context += 1
                context_depths[context] = 1

            if prev == '.':
                dots_end = i-1
                bulges = add_bulge(bulges, (dots_start-1, dots_end+1), context, "4")
            """
            if prev == ')':
                bulges = add_bulge(bulges, (-1,-1), context, "3")
            """

        if brackets[i] == ')':
            if len(opens) == 0:
                error_exit("ERROR: Unmatched close bracket")

            stem_pairs.append((opens.pop(), i))

            context_depths[context] -= 1

            if context_depths[context] == 0:
                finished_bulges += bulges[context]
                bulges[context] = []
                context -= 1
 

            if prev == '.':
                dots_end = i-1
                bulges = add_bulge(bulges, (dots_start-1, dots_end+1), context, "2")
  
        if brackets[i] == '.':
            if prev == '.':
                continue

            dots_start = i

        prev = brackets[i]
    if prev == '.':
        dots_end = i
        bulges = add_bulge(bulges, (dots_start-1, dots_end), context, "7")
    elif prev == '(':
        print >>sys.stderr, "Unmatched bracket at the end"
        sys.exit(1)
    """
    elif prev == ')':
        bulges = add_bulge(bulges, (i+1, i), context, "8")
    """
    
    if context in bulges.keys():
        finished_bulges += bulges[context]

    if len(opens) > 0:
        error_exit("ERROR: Unmatched open bracket")

    stem_pairs.sort()
    stems = condense_stem_pairs(stem_pairs)
    
    return (finished_bulges, stems)

class BulgeGraph:
    """
    A class that stores the information about an RNA graph
    """
    def __init__(self, filename=None):
        self.name = 'unnamed'
        self.seq = ''
        self.length = 0

        # Look up tables
        self.sampled_stems = dict()
        self.coords = dict()
        self.twists = dict()
        self.stem_invs = dict()
        self.bases = dict()
        self.edges = c.defaultdict(set)
        self.defines = dict()

        self.vposs = c.defaultdict( dict )
        self.vvecs = c.defaultdict( dict )
        self.v3dposs = c.defaultdict( dict )
        self.vbases = c.defaultdict( dict )
        self.vinvs = c.defaultdict( dict )

        self.weights = dict()

        # for creating new vertex indeces
        self.vertex_counter = 0
        self.name_counter = 0

        # store which bulges were merged
        self.merge_defs = dict()
        self.longrange = c.defaultdict( set )
        
        self.bp_distances = None

        if filename != None:
            self.parse_graph(filename)

        #self.calc_bp_distances()

    # get an internal index for a named vertex
    # this applies to both stems and edges
    def get_vertex(self, name = None):
        if name == None:
            name = "x%d" % (self.name_counter)
            self.name_counter += 1

        return name

    def has_connection(self, v1, v2):
        """ Is there an edge between these two nodes """

        if v2 in self.edges[v1]:
            return True
        else:
            return False

    def stem_length(self, key):
        d = self.defines[key]

        if key[0] != 's' and len(d) == 4:
            return min(d[1] - d[0], d[3] - d[2]) + 1

        return (d[1] - d[0]) + 1

    def get_node_from_residue_num(self, base_num):
        """
        Iterate over the defines and see which one encompasses this base.
        """
        for key in self.defines.keys():
            define = self.defines[key]

            for i in range(0, len(define), 2):
                a = [int(define[i]), int(define[i+1])]
                a.sort()
                
                if key[0] != 's':
                    if a[0] != 0:
                        a[0] += 1
                    if a[1] != self.length:
                        a[1] -= 1

                if base_num >= a[0] and base_num <= a[1]:
                    return key

        raise Exception("Base number %d not found in the defines." % (base_num))

    def get_stem_angle(self, s1, s2):
        '''
        Return the angle between the vectors of the stems.

        @param s1: The name of the first stem.
        @param s2: The name of the second stem.
        @return: The angle (in radians), between the two stem vectors.
        '''
        s1_vec = self.coords[s1][1] - self.coords[s1][0]
        s2_vec = self.coords[s2][1] - self.coords[s2][0]

        angle = cuv.vec_angle(s1_vec, s2_vec)

        return min(angle, math.pi - angle)
        #return angle

    def get_bulge_angle_stats_core(self, define, connections):
        '''
        Return the angle stats for a particular bulge. These stats describe the
        relative orientation of the two stems that it connects.

        @param define: The name of the bulge.
        @param connections: The two stems that are connected by it.
        @return: cbs.AngleStat object
        '''
        (stem1, twist1, stem2, twist2, bulge) = cgg.get_stem_twist_and_bulge_vecs(self, define, connections)

        (s1b, s1e) = self.get_sides(connections[0], define)
        (s2b, s1e) = self.get_sides(connections[1], define)

        # Get the orientations for orienting these two stems
        (r, u, v, t) = cgg.get_stem_orientation_parameters(stem1, twist1, stem2, twist2)
        (r1, u1, v1) = cgg.get_stem_separation_parameters(stem1, twist1, bulge)

        dims =self.get_bulge_dimensions(define)

        angle_stat = cbs.AngleStat(self.name, dims[0], dims[1], u, v, t, r1, u1, v1, s1b, s2b, self.defines[define])

        return angle_stat

    def virtual_residues(self):
        '''
        Iterate over every virtual residue in the structure.
        '''
        for s in self.stems():
            for i in range(self.stem_length(s)):
                yield((s,i))

    def get_long_range_constraints(self):
        seen = collections.defaultdict(set)
        constr = []

        for key1 in self.longrange.keys():
            for key2 in self.longrange[key1]:
                if key2 in seen[key1]:
                    continue

                seen[key1].add(key2)
                seen[key2].add(key1)

                point1 = self.get_point(key1)
                point2 = self.get_point(key2)

                constr += [(key1, key2, cuv.vec_distance(point1, point2))]

        return constr

    def connections(self, bulge):
        '''
        Return the edges that connect to a bulge in a list form,
        sorted by lowest res number of the connection.
        '''
        connections = list(self.edges[bulge])
        connections.sort(key=lambda x: self.defines[x][0])

        return connections

    def get_bulge_angle_stats(self, bulge):
        '''
        Return the angle stats for a particular bulge. These stats describe the
        relative orientation of the two stems that it connects.

        @param bulge: The name of the bulge.
        @param connections: The two stems that are connected by it.
        @return: The angle statistics in one direction and angle statistics in the other direction
        '''

        if bulge == 'start':
            return (cbs.AngleStat(), cbs.AngleStat())

        #print "bulge:", bulge
        connections = self.connections(bulge)

        angle_stat1 = self.get_bulge_angle_stats_core(bulge, connections)
        angle_stat2 = self.get_bulge_angle_stats_core(bulge, list(reversed(connections)))

        return (angle_stat1, angle_stat2)

    def get_stem_stats(self, stem):
        '''
        Calculate the statistics for a stem and return them. These statistics will describe the
        length of the stem as well as how much it twists.

        @param stem: The name of the stem.

        @return: A StemStat structure containing the above information.
        '''

        ss = cbs.StemStat()

        ss.pdb_name = self.name
        ss.bp_length = abs(self.defines[stem][0] - self.defines[stem][1])
        ss.phys_length = cuv.magnitude(self.coords[stem][0] - self.coords[stem][1])
        ss.twist_angle = cgg.get_twist_angle(self.coords[stem], self.twists[stem])
        ss.define = self.defines[stem]

        return ss

    def are_adjacent_stems(self, s1, s2):
        '''
        Check whether two stems are separated by one intermediate element.

        The intermediate element must be an interior loop.

        For example.

        True: s1-b-s2
        False: s1-b-sx-b-s2

        @param: s1 the name of the first stem
        @param: s2 the name of the second stem
        '''

        for edge in self.edges[s1]:
            # the intermediate element must be an interior loop
            if self.weights[edge] == 2:
                if s2 in self.edges[edge]:
                    return True

        return False

    def are_any_adjacent_stems(self, s1, s2):
        '''
        Check whether two stems are separated by one intermediate element.

        For example.

        True: s1-b-s2
        False: s1-b-sx-b-s2

        @param: s1 the name of the first stem
        @param: s2 the name of the second stem
        '''

        for edge in self.edges[s1]:
            if s2 in self.edges[edge]:
                return True

        return False


    def get_random_bulge(self):
        '''
        Return the name of a random bulge.

        @return: The name of a bulge.
        '''

        bulges = []

        for d in self.defines.keys():
            if d[0] != 's' and len(self.edges[d]) == 2:
                bulges += [d]

        return random.choice(bulges)

    def get_length(self, vertex):
        '''
        Get the minimum length of a vertex.

        If it's a stem, then the result is its length (in base pairs).

        If it's a bulge, then the length is the smaller of it's dimensions.

        @param vertex: The name of the vertex.
        '''
        if vertex[0] == 's':
            return abs(self.defines[vertex][1] - self.defines[vertex][0]) + 1
        else:
            dims = self.get_bulge_dimensions(vertex)

            if dims[0] == 0:
                return dims[1]
            else:
                return dims[0]

    def calc_vres_distance(self, s1, i1, s2, i2):
        '''
        Calculate the number of virtual residues that separate two
        virtual residues.

        @param s1: The first stem
        @param i1: The index into the first stem
        @param s2: The second stem
        @param i2: The index into the second stem
        '''
        side1 = self.closest_sides[s1][s2]
        side2 = self.closest_sides[s2][s1]

        if s1 == s2:
            return abs(i2 - i1)

        if side1 == 0:
            add1 = i1
        else:
            add1 = self.stem_length(s1) - i1

        if side2 == 0:
            add2 = i2
        else:
            add2 = self.stem_length(s2) - i2

        return self.bp_distances[s1][s2] + add1 + add2

    def calc_bp_distances(self):
        '''
        Calculate the least number of bases that separate two elements.

        If they are connected, the result is 0. If they are separated by a
        stem, then the result is the length of the stem. If separated by a bulge
        the distance is min(dims(bulge))

        The values are calculated using a slight variation of the Floyd Warshall Algorithm.

        The result is stored in the double dict self.bp_distances. 
        '''

        dist = cuds.DefaultDict(cuds.DefaultDict(0))
        sides = cuds.DefaultDict(cuds.DefaultDict((0,1)))
        stem_sides = c.defaultdict(lambda: c.defaultdict(lambda: 0))

        defs = self.defines.keys()

        for i in defs:
            for j in defs:
                dist[i][j] = 10000000
                
                if i == j:
                    dist[i][j] = 0

                if j in self.edges[i]:
                    dist[i][j] = 0
                    
                    # stem_sides stores which side of stem i points
                    # to stem j 
                    # if i or j is a bulge, then the stored side is
                    # for the stem
                    if i[0] == 's':
                        stem_sides[i][j] = self.get_sides(i, j)[0]
                        stem_sides[j][i] = self.get_sides(i, j)[0]
                    else:
                        stem_sides[j][i] = self.get_sides(j, i)[0]
                        stem_sides[i][j] = self.get_sides(j, i)[0]

        sys.stderr.write('Calculating base pair distances')
        for k in defs:
            sys.stderr.write('.')
            sys.stderr.flush()
            for i in defs:
                for j in defs:
                    inter_distance = self.get_length(k)

                    if k[0] == 's':
                        # if we want to go through a stem, we have to make sure
                        # that the two nodes being connected are on opposide sides of it
                        if stem_sides[i][k] == stem_sides[j][k]:
                            inter_distance = 0
                    
                    if dist[i][j] > dist[i][k] + inter_distance + dist[k][j]:
                        dist[i][j] = dist[i][k] + inter_distance + dist[k][j]
                        #sides[i][j] = (s1b, s2b)

                        stem_sides[i][j] = stem_sides[i][k]
                        stem_sides[j][i] = stem_sides[j][k]

        sys.stderr.write('Finished\n')
        self.bp_distances = dist
        self.closest_sides = stem_sides

    def breadth_first_traversal(self, start=None):
        '''
        Do a breadth-first traversal of the graph.

        @param start: The starting node.
        @return: The path of nodes visited in the bfs traversal. Returned as a list of triples (stem, bulge, stem).
        '''
        visited = set()
        to_visit = []
        path = []

        if start == None:
            for define in self.defines.keys():
                if define[0] != 's' and self.weights[define] == 1 and len(self.edges[define]) == 1:
                    to_visit.append(('start', 'start', list(self.edges[define])[0]))
        else:
            to_visit.append(('start', 'start', start))

        while len(to_visit) > 0:
            (prev_stem, prev_node, curr_node) = to_visit.pop()

            if curr_node in visited:
                continue

            visited.add(curr_node)

            if curr_node[0] == 's': 
                path += [(prev_stem, prev_node, curr_node)]
                for edge in self.edges[curr_node]:
                    to_visit.append((curr_node, curr_node, edge))
            else:
                if len(self.edges[curr_node]) == 2:
                    for edge in self.edges[curr_node]:
                        to_visit.append((prev_stem, curr_node, edge))

        return path

    def get_bulge_dimensions(self, bulge):
        '''
        Return the dimensions of the bulge.

        If it is single stranded it will be (0, x). Otherwise it will be (x, y).

        @param bulge: The name of the bulge.
        @return: A pair containing its dimensions
        '''

        bd = self.defines[bulge]
        prev_stem = self.connections(bulge)[0]
        c = self.connections(bulge)

        (s1b, s1e) = self.get_sides(c[0], bulge)
        (s2b, s2e) = self.get_sides(c[1], bulge)

        if len(bd) == 2:
            dims = (abs(bd[1] - bd[0]), 0)
        else:
            dims = (abs(bd[1] - bd[0]), abs(bd[3] - bd[2]))

        if s1b == s2b:
            assert(s1b == 0)
            if bd[0] == self.defines[prev_stem][0]:
                return dims
            else:
                assert(bd[0] == self.defines[prev_stem][3])
                return (dims[1], dims[0])
        else:
            if bd[0] == self.defines[prev_stem][1]:
                return dims
            else:
                '''
                cud.pv('bulge')
                cud.pv('(bd[0], bd[1])')
                cud.pv('self.defines[c[0]]')
                cud.pv('self.defines[c[1]]')
                cud.pv('(s1b, s2b)')
                '''
                assert(bd[1] == self.defines[prev_stem][2])
                return (dims[1], dims[0])

    # internal function for creating the forward
    # and reverse edge between two vertices
    def connect_vertices(self, v1, v2):
        if v1 not in self.edges:
            self.edges[v1] = set()

        if v2 not in self.edges:
            self.edges[v2] = set()

        self.edges[v1].add(v2)
        self.edges[v2].add(v1)


    # Add an edge between two named vertices
    def add_edge(self, start, end, weight=1):
        v1 = self.get_vertex(start)
        v2 = self.get_vertex(end)

        self.connect_vertices(v1, v2)


    # Merging is done only two vertices at a time
    # We want to consolidate the definitions so that each merged vertex
    # is defined by all of the vertices that have been consolidated to create it
    # i.e.
    #
    # x0 b1 b2
    # x1 b3 x0
    #
    # becomes
    #
    # x1 b3 b1 b3
    def collapse_merge_defs(self):
        done = False
        merge_defs = self.merge_defs

        while not done:
            done = True
            for key1 in merge_defs.keys():
                for key2 in merge_defs.keys():
                    values = merge_defs[key2]

                    for k in range(len(values)):
                        if values[k] == key1:
                            del values[k]
                            values += merge_defs[key1]

                            del merge_defs[key1] 
                            done = False
    def stems(self):
        for d in self.defines.keys():
            if d[0] == 's':
                yield d

    def elements(self):
        for d in self.defines.keys():
            yield d

    def junctions(self):
        for d in self.defines.keys():
            if d[0] != 's' and self.weights[d] == 1 and len(self.edges[d]) == 2:
                yield d
    
    def stem_like(self):
        for d in self.defines.keys():
            if self.weights[d] == 2 or d[0] == 's':
                yield d
    
    def bulges(self):
        for d in self.defines.keys():
            if len(self.edges[d]) == 2:
                yield d

    def interior_loops(self):
        for d in self.defines.keys():
            if d[0] != 's' and self.weights[d] == 2:
                yield d

    def get_centers(self):
        centers = []
        for d in self.defines.keys():
            if d in self.coords.keys():
                p = self.coords[d]
                centers += [(p[1] + p[0]) / 2.]

        return centers

    def get_name_str(self):
        return "name " + self.name + '\n'

    def get_seq_str(self):
        return "seq " + self.seq + "\n"

    def get_length_str(self):
        return "length " + str(self.length) + '\n'

    def print_name(self):
        print self.get_name_str()


    def get_point(self, key):
        coords = self.coords[key]

        if key[0] == 'b':
            if len(self.edges[key]) == 1:
                # the loops centroid is always the second point
                return coords[1]
            else:
                # otherwise take the middle of the bulge
                return (coords[0] + coords[1]) / 2.0
        else:
            # we'll consider the middle of the stem
            return (coords[0] + coords[1]) / 2.0


    def get_sides_plus(self, s1, b):
        '''
        Get the side of s1 that is next to b.

        s1e -> s1b -> b

        @param s1: The stem.
        @param b: The bulge.
        @return: A tuple indicating which side is the one next to the bulge
                 and which is away from the bulge.
        '''
        s1d = self.defines[s1]
        bd = self.defines[b]

        #print >>sys.stderr, "s1: %s b: %s" % (s1, b)

        for i in xrange(4):
            for k in xrange(len(bd)):
                if s1d[i] == bd[k]:
                    return (i, k)
        return None

    def get_sides(self, s1, b):
        '''
        Get the side of s1 that is next to b.

        s1e -> s1b -> b

        @param s1: The stem.
        @param b: The bulge.
        @return: A tuple indicating which side is the one next to the bulge
                 and which is away from the bulge.
        '''
        s1d = self.defines[s1]
        bd = self.defines[b]

        #print >>sys.stderr, "s1: %s b: %s" % (s1, b)

        for i in xrange(4):
            for k in xrange(len(bd)):
                if s1d[i] == bd[k]:
                    if i == 0 or i == 3:
                        s1b = 0
                    else:
                        s1b = 1

        if s1b == 0:
            s1e = 1
        else:
            s1e = 0

        return (s1b, s1e)

    def get_define_str(self):
        defines_str = ''
        for key in self.defines.keys():
            defines_str += "define %s %d %s" % ( key, self.weights[key], " ".join([str(d) for d in self.defines[key]]))
            defines_str += '\n'
        return defines_str

    def print_defines(self):
        print self.get_define_str()

    # Print the definitions for the merged nodes
    def print_merge_defs(self):
        self.collapse_merge_defs()

        for key in self.merge_defs.keys():
            out_str = "merge %s" % (key)

            for val in self.merge_defs[key]:
                out_str += " %s" % (key)
            print out_str
    
    # Take the define lines, combine them with the merge definitions
    # and output the actual merged coords of each motif
    def print_merged_defines(self, defines):
        self.collapse_merge_defs()
        merged = set()

        for key in self.merge_defs.keys():
            out_str = "define %s" % (key)

            for val in self.merge_defs[key]:
                name = val
                merged.add(name)
                orig_define = " ".join(defines[name])
                out_str += " %s" % (orig_define)

            print out_str
        
        for define in defines.keys():
            if define not in merged:
                print "define %s %s" % (define, " ".join(defines[define]))
            
    def get_connect_str(self):
        whole_str = ''
        for key in self.edges:
            # Our graph will be defined by the stems and the bulges they connect to
            name = key
            if name[0] == 's':
                out_str = "connect %s" % (name)

                for dest in self.edges[key]:
                    out_str += " %s" % (dest)

                whole_str += out_str
                whole_str += '\n'
        return whole_str


    def rotate_coords(self, rotation):
        for key in self.coords.keys():
            coords = self.coords[key]

            new_coords = (np.dot(rotation, coords[0]), np.dot(rotation, coords[1]))

            self.coords[key] = new_coords
            if key[0] == 's':
                twists = self.twists[key]
                new_twists = (np.dot(rotation, twists[0]), np.dot(rotation, twists[1]))
                self.twists[key] = new_twists

    def get_twists(self, node):
        '''
        Get the array of twists for this node. If the node is a stem,
        then the twists will simply those stored in the array.

        If the node is an interior loop or a junction segment, 
        then the twists will be the ones that are adjacent to it.

        If the node is a hairpin loop or a free end, then only
        one twist will be returned.

        @param node: The name of the node
        '''
        if node[0] == 's':
            return self.twists[node]

        connections = list(self.edges[node])
        (s1b, s1e) = self.get_sides(connections[0], node)

        if len(connections) == 1:
            return (cuv.normalize(cuv.vector_rejection(
                             self.twists[connections[0]][s1b],
                             self.coords[connections[0]][1] - 
                             self.coords[connections[0]][0])),)

        if len(connections) == 2:
            # interior loop or junction segment
            (s2b, s2e) = self.get_sides(connections[1], node)

            bulge_vec = (self.coords[connections[0]][s1b] - 
                         self.coords[connections[1]][s2b]) 

            return (cuv.normalize(cuv.vector_rejection(
                    self.twists[connections[0]][s1b],
                    bulge_vec)),
                    cuv.normalize(cuv.vector_rejection(
                    self.twists[connections[1]][s2b],
                    bulge_vec)))

        # uh oh, this shouldn't happen since every node
        # should have either one or two edges
        return None

    def translate_coords(self, translation):
        '''
        Translate all of the coordinates a certain distance.

        The twists are always relative to the stems, so they do not need to be
        translated.

        @param translation. The translation vector.
        '''
        for key in self.coords.keys():
            coords = self.coords[key]

            new_coords = (coords[0] + translation, coords[1] + translation)

            self.coords[key] = new_coords

    # Print the graph in the form
    # 
    # s1 b1 b2
    # 
    # Which would imply that s1 is connected to b1 and b2 
    def print_me(self):
        print self.get_connect_str()

    def get_sampled_stems_str(self):
        out_str = ''
        for key in self.sampled_stems.keys():
            out_str += 'sampled %s %s\n' % (key, " ".join(map(str, self.sampled_stems[key])))
        return out_str

    def get_coord_str(self):
        out_str = ''
        for key in self.coords.keys():
            [p, n] = self.coords[key]
            out_str += "coord %s %s %s" % (key, " ".join([str(pt) for pt in p]), " ".join([str(pt) for pt in n]))
            out_str += '\n'
        return out_str
    
    def get_twist_str(self):
        out_str = ''
        for key in self.twists.keys():
            [p, n] = self.twists[key]
            out_str += "twist %s %s %s" % (key, " ".join([str(pt) for pt in p]), " ".join([str(pt) for pt in n]))
            out_str += '\n'
        return out_str

    def get_long_range_str(self):
        out_str = ''
        for key1 in self.longrange.keys():
            for key2 in self.longrange[key1]:
                out_str += "longrange %s %s\n" % (key1, key2)

        return out_str

    def print_coords(self):
        print self.get_coord_str()

    def print_twists(self):
        print self.get_twist_str()

    def get_as_str(self):
        out_str = ''
        out_str += self.get_name_str()
        out_str += self.get_seq_str()
        out_str += self.get_length_str()
        out_str += self.get_define_str()
        out_str += self.get_coord_str()
        out_str += self.get_sampled_stems_str()
        out_str += self.get_twist_str()
        out_str += self.get_connect_str()
        out_str += self.get_long_range_str()

        return out_str
    
    def copy(self):
        bg = BulgeGraph()

        bg.name = self.name
        bg.seq = self.seq
        bg.length = self.length

        bg.sampled_stems = self.sampled_stems.copy()
        bg.coords = self.coords.copy()
        bg.twists = self.twists.copy()
        bg.edges = self.edges.copy()
        bg.defines = self.defines.copy()
        bg.weights = self.weights.copy()
        bg.merge_defs = self.merge_defs.copy()
        bg.longrange = self.longrange.copy()

        if self.bp_distances is not None:
            bg.bp_distances = self.bp_distances.copy()

        return bg

    def output(self, filename):
        f = open(filename, 'w')
        out_str = self.get_as_str()

        f.write(out_str)

        f.close()

    def create_bulge_graph(self, stems, bulges):
        '''
        Find out which stems connect to which bulges

        Stems and bulges which share a nucleotide are considered connected.

        :param stems: A list of tuples of tuples of the form [((s1, e1), (s2, e2))]
                      where s1 and e1 are the nucleotides at one end of the stem
                      and s2 and e2 are the nucleotides at the other.

        :param bulges: A list of tuples of the form [(s, e)] where s and e are the 
                       numbers of the nucleotides at the start and end of the bulge.
        '''
        for i in range(len(stems)):
            stem = stems[i]
            for j in range(len(bulges)):
                bulge = bulges[j]
                if any_difference_of_one(stem, bulge):
                    self.edges['s%d' % (i)].add('b%d' % (j))
                    self.edges['b%d' % (j)].add('s%d' % (i))

    def get_flanking_region(self, bulge_name, side=0):
        '''
        If a bulge is flanked by stems, return the lowest residue number
        of the previous stem and the highest residue number of the next
        stem.

        @param bulge_name: The name of the bulge
        @param side: The side of the bulge (indicating the strand)
        '''
        d = (self.defines[bulge_name][2*side], self.defines[bulge_name][2*side+1])
        m1 = d[0]
        m2 = d[1]

        for edge in self.edges[bulge_name]:
            # should always be a stem...
            d1 = self.defines[edge]

            if d[0] == d1[1]:
                m1 = d1[0]

            if d[0] == d1[3]:
                m1 = d1[2]

            if d[1] == d1[0]:
                m2 = d1[1]

            if d[1] == d1[2]:
                m2 = d1[3]

        if m1 == 0:
            m1 = 1

        return (m1, m2)
        

    def get_flanking_sequence(self, bulge_name, side=0):
        if len(self.seq) == 0:
            print >>sys.stderr, "No sequence present in the bulge_graph: %s" % (self.name)

        (m1, m2) = self.get_flanking_region(bulge_name, side)

        return self.seq[m1-1:m2]

    def get_flanking_handles(self, bulge_name, side=0):
        '''
        Get the indices of the residues for fitting bulge regions.
        
        So if there is a loop like so (between residues 7 and 16):

        (((...))))
        7890123456
          ^   ^  

        Then residues 9 and 13 will be used as the handles against which
        to align the fitted region.

        In the fitted region, the residues (2,6) will be the ones that will
        be aligned to the handles.

        @return: (orig_chain_res1, orig_chain_res1, flanking_res1, flanking_res2)
        '''
        def1 = self.defines[bulge_name]
        f1 = self.get_flanking_region(bulge_name, side)

        a = def1[side*2]
        b = def1[side*2+1]

        i1 = def1[side*2] - f1[0]
        i2 = def1[side*2 + 1] - f1[0]

        #cud.pv('i1')
        #cud.pv('i2')
        '''
        if b == self.length:
            b -= 1
            i2 -= 1

        #return (def1[side*2], def1[side*2+1], def1[side*2] - f1[0], def1[side*2 + 1] - f1[0])
        '''
        return (a, b, i1, i2)


    def dump(self):
        print self.get_as_str()

    def create_stem_graph(self, stems, bulge_counter):
        '''
        Determine which stems are connected to each other. A stem can be connected to
        another stem when there is an interior loop with an unpaired nucleotide on
        one side. In this case, a bulge will be created on the other side, but it
        will only consist of the two paired bases around where the unpaired base 
        would be if it existed.

        The defines for these bulges will be printed as well as the connection strings
        for the stems they are connected to.

        :param stems: A list of tuples of tuples of the form [((s1, e1), (s2, e2))]
                      where s1 and e1 are the nucleotides at one end of the stem
                      and s2 and e2 are the nucleotides at the other.
        :param bulge_counter: The number of bulges that have been encountered so far.

        :returns: A dictionary indexed by the number of a stem, containing a set of the 
                 other stems that the index is connected to.
        '''
        #print "stems:", stems
        stem_stems = dict()
        define_text = ""
        connect_text = ""
        for i in range(len(stems)):
            stem = stems[i]
            for j in range(i+1, len(stems)):
                for k1 in range(2):
                    # don't fear the for loop
                    for k2 in range(2):
                        for l1 in range(2):
                            for l2 in range(2):
                                s1 = stems[i][k1][l1]
                                s2 = stems[j][k2][l2]
                                if abs(s1 - s2) == 1:
                                    stem_stems_set = stem_stems.get(i, set())
                                    if j not in stem_stems_set:
                                        bn = 'b%d' % (bulge_counter)
                                        self.defines[bn] = [min(s1, s2)+1, max(s1, s2)+1]
                                        self.weights[bn] = 1

                                        self.edges['s%d' % i].add(bn)
                                        self.edges[bn].add('s%d' % i)

                                        self.edges['s%d' % j].add(bn)
                                        self.edges[bn].add('s%d' % j)

                                        bulge_counter += 1
                                        stem_stems_set.add(j)
                                    stem_stems[i] = stem_stems_set
        return stem_stems

    # Delete a vertex after merging it with another
    def remove_vertex(self, v):
        '''
        Delete a vertex after merging it with another
        '''
        # delete all edges to this node
        for key in self.edges[v]:
            self.edges[key].remove(v)

        for edge in self.edges:
            if v in self.edges[edge]:
                self.edges[edge].remove(v)

        # delete all edges from this node
        del self.edges[v]
        del self.defines[v]

    def reduce_defines(self):
        """
        Make defines like this:

        define x0 2 124 124 3 4 125 127 5 5
        
        Into this:

        define x0 2 3 5 124 127

        That is, consolidate contiguous bulge region defines.
        """

        new_defines = []

        for key in self.defines.keys():
            if key[0] != 's':
                assert(len(self.defines[key]) % 2 == 0)

                j = 0
                new_j = 0
                #print >>sys.stderr,"self.defines[key]:", self.defines[key]

                while new_j < len(self.defines[key]):

                    j = new_j
                    new_j += j + 2

                    #print >>sys.stderr, "reduce_defines", i, j
                    (f1, t1) = (int(self.defines[key][j]), int(self.defines[key][j+1]))

                    # remove bulges of length 0
                    if f1 == -1 and t1 == -2:
                        del self.defines[key][j]
                        del self.defines[key][j]
                        
                        new_j = 0
                        continue

                    # merge contiguous bulge regions
                    for k in range(j+2, len(self.defines[key]), 2):
                        (f2, t2) = (int(self.defines[key][k]), int(self.defines[key][k+1]))

                        #print >>sys.stderr, "j: %d f1: %d, t1: %d, f2: %d, t2: %d" % (j, f1, t1, f2, t2)

                        if t2 + 1 != f1 and t1 + 1 != f2:
                            continue

                        #print >>sys.stderr, "pre self.defines[key]:", self.defines[key]

                        if t2 + 1 == f1:
                            self.defines[key][j] = str(f2)
                            self.defines[key][j+1] = str(t1)
                        elif t1 + 1 == f2:
                            self.defines[key][j] = str(f1)
                            self.defines[key][j+1] = str(t2)

                        del self.defines[key][k]
                        del self.defines[key][k]

                        #print >>sys.stderr, "post self.defines[key]:", self.defines[key]
                        
                        new_j = 0

                        break

    def merge_vertices(self, vertices):
        '''
        This is done when two of the outgoing strands of a stem
        go to different bulges
        It is assumed that the two ends are on the same sides because
        at least one vertex has a weight of 2, implying that it accounts
        for all of the edges going out of one side of the stem

        :param vertices: A list of vertex names to combine into one.
        '''
        merge_str = ""
        new_vertex = self.get_vertex()
        self.weights[new_vertex] = 0

        #assert(len(vertices) == 2)

        connections = set()
        needs_merging = set()

        for v in vertices:
            merge_str += " %s" % (v)

            # what are we gonna merge?
            for item in self.edges[v]:
                connections.add(item)

            # Add the definition of this vertex to the new vertex
            #self.merge_defs[new_vertex] = self.merge_defs.get(new_vertex, []) + [v]

            if v[0] == 's':
                self.defines[new_vertex] = self.defines.get(new_vertex, []) + [self.defines[v][0], self.defines[v][2]] + [self.defines[v][1], self.defines[v][3]]
            else:
                self.defines[new_vertex] = self.defines.get(new_vertex,[]) + self.defines[v]


            self.weights[new_vertex] += 1

            # remove the old vertex, since it's been replaced by new_vertex
            self.remove_vertex(v)
            self.reduce_defines()

        #self.weights[new_vertex] = 2

        for connection in connections:
            self.connect_vertices(new_vertex, connection)

        return new_vertex

    def find_bulge_loop(self, vertex):
        '''
        Find a set of nodes that form a loop containing the
        given vertex and being no greater than 4 nodes long.

        :param vertex: The vertex to start the search from.
        :returns: A list of the nodes in the loop.
        '''
        visited = set()
        to_visit = [(key, 0) for key in self.edges[vertex]]
        visited.add(vertex)
        in_path = [vertex]

        while len(to_visit) > 0:
            (current, depth) = to_visit.pop()
            visited.add(current)

            in_path = in_path[:depth]
            in_path.append(current)

            for key in self.edges[current]:
                if key == vertex and depth > 1:
                    if len(in_path[:depth+1]) > 4:
                        continue
                    else:
                        return in_path[:depth+1]
                
                if key not in visited:
                    to_visit.append((key, depth+1))
        return []

    def collapse(self):
        '''
        If any vertices form a loop, then they are either a bulge region of a fork region. The bulge (interior loop) regions will be condensed into one node.
        '''
        # for all stems of length 1, merge with adjacent bulges
        for key in self.edges.keys():
            if key[0] == 's':
                if self.stem_length(key) == 1:
                    self.dissolve_stem(key)

        for key in self.edges.keys():
            if key[0] == 's':
                while True:
                    loop = self.find_bulge_loop(key)
                    bulge_loop = [v for v in loop if v[0] == 'b' or v[0] == 'x']

                    if len(bulge_loop) > 0:
                        #assert(len(bulge_loop) != 1)
                        self.merge_vertices(bulge_loop)

                    if len(bulge_loop) == 0:
                        break

    def print_edges_by_name(self, v):
        '''
        Debug function for printing the connections of a particular vertex
        '''
        output_str = ""
        for item in self.edges[v]:
            output_str += " || %s - %s " % (name, item[0])
        print output_str

    def get_bulge_vertices(self):
        """
        Get all of the vertices which define a bulge. 

        A bulge vertex can only have two connections: to the two stems which it links. 
        Futhermore, each of these stems can
        """

        vertices = []

        for vertex in self.edges:
            if self.weights[vertex] == 2:
                vertices += [vertex]
        
        return vertices

    def get_bulge_vertex_names(self):
        bulges = self.get_bulge_vertices()

        return [v for v in bulges]

    def get_named_define(self, vertex_name):
        #print "vertex_name:", vertex_name
        #print "defines:", self.defines

        return self.defines[vertex_name]

    def get_bulged_stems(self):
        '''
        Get a list of the stem pairs which are separated by a bulge
        region.
        '''
        bulges = self.get_bulge_vertices()
        stem_pairs = []

        for bulge in bulges:
            connections = self.edges[bulge]
            stem_pairs += [tuple(connections)]

        return stem_pairs

    def get_bulged_stem_names(self):
        stem_pairs = self.get_bulged_stems()
        stem_pair_names = []

        for stem_pair in stem_pairs:
            stem_pair_names += [(stem_pair[0], stem_pair[1])]

        return stem_pair_names

    def add_node(self, name, edges, define, weight=1):
        self.defines[name] = define
        self.edges[name] = edges
        self.weights[name] = weight

        for edge in self.edges[name]:
            self.edges[edge].add(name)


    def dissolve_stem(self, key):
        #print >>sys.stderr, "hi"
        d = self.defines[key]

        bulge_sides = dict()
        bulge_sides[0] = set()
        bulge_sides[1] = set()

        #print >>sys.stderr,"dissolving:", key
        #print >>sys.stderr, "edges:", self.edges[key]
        for edge in self.edges[key]:
            be = self.defines[edge]

            for i in range(0, 4):
                for j in range(0, len(be)):
                    if d[i] == be[j]:
                        #print >>sys.stderr, key, edge
                        bulge_sides[i / 2].add(edge)
        
        #print >>sys.stderr, "bulge_sides:", bulge_sides
        
        new_nodes = [0,0]
        for i in range(2):
            new_node = self.get_vertex()
            edges = set()

            mins = 10000
            maxs = -10000

            for bulge in bulge_sides[i]:
                if bulge in self.defines:
                    bd = self.defines[bulge]
                else:
                    continue
                
                for edge in self.edges[bulge]:
                    edges.add(edge)

                mins = min(mins, min(bd))
                maxs = max(maxs, max(bd))
                self.remove_vertex(bulge)

            edges.remove(key)
            #print >> sys.stderr, "new_node", new_node, "edges:", edges, "mins:", mins, "maxs:", maxs
           
            self.add_node(new_node, edges, [mins, maxs], self.weights[bulge])
            new_nodes[i] = new_node

        if len(self.edges[new_nodes[0]]) == 1 and len(self.edges[new_nodes[1]]) == 1:
            dnn0 = self.defines[new_nodes[0]]
            dnn1 = self.defines[new_nodes[1]]
            newer_node = self.get_vertex()

            define = [min(dnn0[0], dnn1[0]), max(dnn0[1], dnn1[1])]
            edges = self.edges[new_nodes[0]].union(self.edges[new_nodes[1]])

            self.add_node(newer_node, edges, define, self.weights[new_nodes[0]])

            self.remove_vertex(new_nodes[0])
            self.remove_vertex(new_nodes[1])

        self.remove_vertex(key)


    # If any vertices form a loop, then they are either a bulge region of a fork region
    def collapse(self):

        # for all stems of length 1, merge with adjacent bulges
        for key in self.edges.keys():
            if key[0] == 's':
                if self.stem_length(key) == 1:
                    self.dissolve_stem(key)

                """
                s1s = int(self.defines[key][0])
                s1e = int(self.defines[key][2])

                if s1s == s1e:
                    #self.dissolve_stem(key)
                    to_merge = [key]

                    bulge = True
                    for key1 in self.edges[key]:
                        if self.weights[key1] != 2:
                            bulge = False
                            #break


                    to_merge += self.edges[key]
                    #print "to_merge:", to_merge
                    new_vertex = self.merge_vertices(to_merge)
                    if bulge:
                        self.weights[new_vertex] = 2
                """


        for key in self.edges.keys():
            if key[0] == 's':
                while True:
                    loop = self.find_bulge_loop(key)
                    bulge_loop = [v for v in loop if v[0] == 'b' or v[0] == 'x']

                    if len(bulge_loop) > 0:
                        #assert(len(bulge_loop) != 1)
                        self.merge_vertices(bulge_loop)

                    if len(bulge_loop) == 0:
                        break

    def parse_graph(self, filename):
        """
        Parse a file and extract the information about the graph.

        @param filename: The name of the file containing the graph information.
        """
        if filename == '-':
            f = open(sys.stdin, 'r')
        else:
            f = open(filename, 'r')

        for line in f:
            if line.strip().find('name') == 0:
                parts = line.strip().split(' ')
                self.name = parts[1]

            if line.strip().find('seq') == 0:
                parts = line.strip().split(' ')
                if len(parts) == 1:
                    self.seq = ''
                else:
                    self.seq = parts[1]

            if line.strip().find('sampled') == 0:
                parts = line.strip().split(' ')
                self.sampled_stems[parts[1]] = [parts[2], int(parts[3]), int(parts[4]), int(parts[5]), int(parts[6])]

            if line.strip().find('coord') == 0:
                parts = line.strip().split(' ')
                self.coords[parts[1]] = [np.array([float(parts[2]), float(parts[3]), float(parts[4])]), np.array([float(parts[5]), float(parts[6]), float(parts[7])])]

            if line.strip().find('twist') == 0:
                parts = line.strip().split(' ')
                self.twists[parts[1]] = [np.array([float(parts[2]), float(parts[3]), float(parts[4])]), np.array([float(parts[5]), float(parts[6]), float(parts[7])])]

            if line.strip().find('length') == 0:
                parts = line.strip().split(' ')
                self.length = int(parts[1])

            if line.strip().find('longrange') == 0:
                parts = line.strip().split(' ')
                for i in range(2, len(parts)):
                    self.longrange[parts[1]].add(parts[i])

            if line.strip().find('connect') == 0:
                parts = line.strip().split(' ')
                source = parts[1]
                stem_num = int(source[1:])
                for i in range(2, len(parts)):
                    self.add_edge(parts[1], parts[i])

            if line.strip().find('define') == 0:
                parts = line.strip().split(' ')
                element = parts[1]

                weight = int(parts[2])
                definition = parts[3:]

                if element in self.defines.keys():
                    print >>sys.stderr, "Already defined:", element
                    sys.exit(1)
                else:
                    vertex = self.get_vertex(element)
                    self.weights[vertex] = weight
                    self.defines[vertex] = [int(d) for d in definition]

                    #print >>sys.stderr, "vertex:", vertex, "element:", element, "definiton:", definition
        f.close()


    def from_dotbracket(self, dotbracket_str):
        (bulges, stems) = find_bulges_and_stems(dotbracket_str)

        for i in range(len(stems)):
            # one is added to each coordinate to make up for the fact that residues are 1-based
            ss1 = stems[i][0][0]+1
            ss2 = stems[i][0][1]+1
            se1 = stems[i][1][0]+1
            se2 = stems[i][1][1]+1

            self.defines['s%d' % (i)] = [min(ss1,se1), max(ss1,se1), 
                                         min(ss2,se2), max(ss2,se2)]
            self.weights['s%d' % (i)] = 1
        for i in range(len(bulges)):
            bulge = bulges[i]
            self.defines['b%d' % (i)] = [bulge[0]+1, bulge[1]+1]
            self.weights['b%d' % (i)] = 1

        self.dotbracket_str = dotbracket_str
        self.seq_length = len(dotbracket_str)
        self.create_bulge_graph(stems, bulges)
        self.create_stem_graph(stems, len(bulges))
        self.collapse()

    def from_dotbracket_file(self, fn):
        f = open(fn, 'r')
        brackets = "".join(f.readlines()).replace('\n', '')
        self.length = len(brackets)
        self.from_dotbracket(brackets)
