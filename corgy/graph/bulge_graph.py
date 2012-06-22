#!/usr/bin/python

import sys
from numpy import array, dot
from data_structures import DefaultDict

def error_exit(message):
    print >>sys.stderr, message
    sys.exit(1)

class BulgeGraph:
    """
    A class that stores the information about an RNA graph
    """
    def __init__(self, filename=None):
        self.name = 'unnamed'

        # Look up tables
        self.coords = dict()
        self.twists = dict()
        self.edges = dict()
        self.defines = dict()
        self.weights = dict()

        # for creating new vertex indeces
        self.vertex_counter = 0
        self.name_counter = 0

        # store which bulges were merged
        self.merge_defs = dict()
        self.length = 0
        self.longrange = DefaultDict(set())

        if filename != None:
            self.parse_graph(filename)

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

    def get_stem_length(self, key):
        d = self.defines[key]
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

    def get_centers(self):
        centers = []
        for d in self.defines.keys():
            if d in self.coords.keys():
                p = self.coords[d]
                centers += [(p[1] + p[0]) / 2.]

        return centers
                        

    def get_name_str(self):
        return "name " + self.name + '\n'

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


    def get_sides(self, s1, b):
        s1d = self.defines[s1]
        bd = self.defines[b]

        #print >>sys.stderr, "s1: %s b: %s" % (s1, b)

        for i in range(4):
            for k in range(len(bd)):
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

            new_coords = (dot(rotation, coords[0]), dot(rotation, coords[1]))

            self.coords[key] = new_coords
            if key[0] == 's':
                twists = self.twists[key]
                new_twists = (dot(rotation, twists[0]), dot(rotation, twists[1]))
                self.twists[key] = new_twists

    def translate_coords(self, translation):
        for key in self.coords.keys():
            coords = self.coords[key]

            new_coords = (coords[0] + translation, coords[1] + translation)

            self.coords[key] = new_coords

            if key[0] == 's':
                twists = self.twists[key]
                new_twists = (twists[0] + translation, twists[1] + translation)
                self.twists[key] = new_twists

    # Print the graph in the form
    # 
    # s1 b1 b2
    # 
    # Which would imply that s1 is connected to b1 and b2 
    def print_me(self):
        print self.get_connect_str()

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
        out_str += self.get_length_str()
        out_str += self.get_define_str()
        out_str += self.get_coord_str()
        out_str += self.get_twist_str()
        out_str += self.get_connect_str()
        out_str += self.get_long_range_str()

        return out_str

    def output(self, filename):
        f = open(filename, 'w')
        out_str = self.get_as_str()

        f.write(out_str)

        f.close()

    def dump(self):
        print self.get_as_str()
        
    # Delete a vertex after merging it with another
    def remove_vertex(self, v):
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


    # This is done when two of the outgoing strands of a stem
    # go to different bulges
    # It is assumed that the two ends are on the same sides because
    # at least one vertex has a weight of 2, implying that it accounts
    # for all of the edges going out of one side of the stem
    def merge_vertices(self, vertices):
        merge_str = ""
        new_vertex = self.get_vertex()
        self.weights[new_vertex] = 0

        #print >>sys.stderr, "Merging:", vertices

        #assert(len(vertices) == 2)

        connections = set()
        needs_merging = set()
        #print >>sys.stderr, "merging:", vertices

        for v in vertices:
            merge_str += " %s" % (v)

            # what are we gonna merge?
            for item in self.edges[v]:
                connections.add(item)

            # Add the definition of this vertex to the new vertex
            #self.merge_defs[new_vertex] = self.merge_defs.get(new_vertex, []) + [v]

            #print >>sys.stderr, "connection_names:", " ".join(connection_names)

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

    # Debug function for printing the connections of a particular vertex
    def print_edges_by_name(self, v):
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
                if self.get_stem_length(key) == 1:
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

            if line.strip().find('coord') == 0:
                parts = line.strip().split(' ')
                self.coords[parts[1]] = [array([float(parts[2]), float(parts[3]), float(parts[4])]), array([float(parts[5]), float(parts[6]), float(parts[7])])]

            if line.strip().find('twist') == 0:
                parts = line.strip().split(' ')
                self.twists[parts[1]] = [array([float(parts[2]), float(parts[3]), float(parts[4])]), array([float(parts[5]), float(parts[6]), float(parts[7])])]

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


