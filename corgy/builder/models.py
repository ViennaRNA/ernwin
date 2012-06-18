#!/usr/bin/python

from corgy.graph.bulge_graph import BulgeGraph
from corgy.utilities.vector import normalize, magnitude, rotation_matrix, vec_angle, get_non_colinear_unit_vector
from corgy.utilities.numeric import get_random_in_range
from corgy.utilities.data_structures import DefaultDict
from corgy.builder.rmsd import centered_rmsd, rmsd, optimal_superposition
from corgy.builder.loops import get_random_spherical
from corgy.visual.pymol import PymolPrinter

from random import random, choice
from math import pi, sin, cos, acos
from numpy import array, cross, dot


avg_stem_bp_length = 2.24

loop_lengths = [ 
        (0., 0.),
        ( 7.0 , 9.0 ), 
        ( 7.459 , 9.33 ), 
        ( 7.774 , 8.945 ), 
        ( 8.102 , 8.985 ), 
        ( 6.771 , 8.182 ), 
        ( 6.465 , 7.533 ), 
        ( 6.435 , 7.676 ), 
        ( 6.605 , 8.987 ), 
        ( 8.396 , 9.367 ), 
        ( 12.13 , 18.68 ), 
        ( 19.76 , 22.32 ), 
        ( 11.57 , 14.59 ), 
        ( 8.702 , 8.744 ), 
        ( 15.46 , 15.46 ), 
        ( 15.0 , 30.0 ), 
        ( 15.0 , 30.0 ), 
        ( 15.0 , 30.0 ), 
        ( 15. , 30. ), 
        ( 15. , 30. ), 
        ( 15. , 30. ), 
        ( 15. , 30. ), 
        ( 15. , 30. ), 
        ( 15. , 30. ), 
        ( 33.02 , 33.02 ) ]


def add_bulge(length, stem1, angle1):
    start = stem1[1]
    stem1_vec = stem1[1] - stem1[0]

    stem1_ncl = get_non_colinear_unit_vector(stem1_vec)
    comp1 = cross(stem1_vec, stem1_ncl)
    comp2 = cross(stem1_vec, comp1)

    rot_mat1 = rotation_matrix(comp2, angle1)
    rot_mat2 = rotation_matrix(stem1_vec, random() * 2.0 * pi)

    new_vec1 = dot(rot_mat1, stem1_vec)
    new_vec2 = dot(rot_mat2, new_vec1)

    nv2 = normalize(new_vec2)
    return (start, start + length * nv2)

def add_stem2(length, stem1, bulge, angle1, angle2):
    stem1_vec = stem1[1] - stem1[0]
    bulge_vec = bulge[1] - bulge[0]
    start = bulge[1]
    
    s1b_cross = cross(bulge_vec, stem1_vec)

    #comp1 = cross(s1b_cross, stem1_vec)
    #comp2 = cross(comp1, stem1_vec)

    #print >>sys.stderr, "angle1: %f angle2: %f" % (angle1, angle2)
    rot_mat1 = rotation_matrix(s1b_cross, angle1)
    #rot_mat2 = rotation_matrix(bulge_vec, -angle2)
    rot_mat2 = rotation_matrix(bulge_vec, pi - angle2)

    new_vec1 = dot(rot_mat1, bulge_vec)
    new_vec2 = dot(rot_mat2, new_vec1)
    #new_vec2 = new_vec1
   
    nv2 = normalize(new_vec2)
    return (start, start + length * nv2)

def get_bulge_dimensions(bg, key):
    bd = bg.defines[key]
    if len(bd) == 2:
        return (0, abs(bd[0] - bd[1]) + 1)
    else:
        dim1 = abs(bd[0] - bd[1]) + 1
        dim2 = abs(bd[2] - bd[3]) + 1
        return (min(dim1, dim2), max(dim1, dim2))

def get_bulge_length(bg, key, length_stats):
    dimensions = get_bulge_dimensions(bg, key)
    bulge_length = dimensions[0]

    #print >>sys.stderr, "dimensions:", dimensions
    found = False
    
    dim1 = dimensions[1]
    dim0 = dimensions[0]
    while not found:
        try:
            length = choice(length_stats['bulge_mid'][dim0][dim1])
            found = True
        except IndexError:
            found = False
            dim1 -= 1

    return length


def get_loop_length(bg, key):
    if int(bg.defines[key][0]) > int(bg.defines[key][1]):
        loop_length = 1. 
    else:
        loop_length = int(bg.defines[key][1]) - int(bg.defines[key][0])
    #print >>sys.stderr, "loop_length:", loop_length

    return loop_lengths[loop_length][0] + random() * (loop_lengths[loop_length][1] - loop_lengths[loop_length][0])


def add_close_segment(length, prev_stem):
    start = prev_stem[1]
    direction = prev_stem[1] - prev_stem[0]

    #print >>sys.stderr, "direction:", direction
    comp1 = get_non_colinear_unit_vector(direction)
    #print >>sys.stderr, "comp1:", comp1
    comp2 = cross(comp1, direction)

    rm1 = rotation_matrix(comp2, random() * pi / 2)
    rm2 = rotation_matrix(direction, random() * 2 * pi)

    nv1 = dot(rm1, direction)
    nv2 = dot(rm2, nv1)

    nv = normalize(nv2)

    return (start, start + length * nv)

class SpatialModel:
    def __init__(self, bg, length_stats, angle_stats):
        self.length_stats = length_stats
        self.angle_stats = angle_stats
        self.bg = bg
        self.segments = DefaultDict(DefaultDict(0))
        self.pymol_printer = PymolPrinter()

        pass


    def find_start_node(self):
        #print >>sys.stderr, "bg.edges", bg.edges
        for define in self.bg.defines.keys():
            if define[0] != 's' and self.bg.weights[define] == 1 and len(self.bg.edges[define]) == 1:
                self.to_visit += [(define, 'start', 1)]
                break

    def get_transform(self, edge):
        segment = self.vectors[edge]
        return segment[0]

    def get_rotation(self, edge):
        segment = self.vectors[edge]

        desired_orientation = array([1., 0., 0.])
        vec = normalize(segment[1] - segment[0])
        norm = cross(desired_orientation, vec)

        angle = vec_angle(desired_orientation, vec)
        return rotation_matrix(norm, angle)

    def add_stem_segment(self, edge, other_edge, segment):
        self.bg.coords[edge] = segment
        """
        if edge == self.rot_segment:
            self.pymol_printer.add_segment(segment[0], segment[1], 'purple', 5.0, edge)
        else:
            self.pymol_printer.add_segment(segment[0], segment[1], 'green', 3.0, edge)
        """

        if self.first_segment:
            self.first_translation = self.get_transform(edge)
            self.first_rotation = self.get_rotation(edge)

            self.first_segment = False

    def add_loop_segment(self, edge, other_edge, segment):
        self.bg.coords[edge] = segment

    def add_interior_bulge(self, edge, other_edge, segment):
        self.bg.coords[edge] = segment

    def add_single_bulge(self, edge, other_edge, segment):
        self.bg.coords[edge] = segment

    def add_loop(self, prev_node, curr_node, side):
            ## Loop
            ## Must be connected to a stem
            ## Add a vector for this loop
            length = get_loop_length(self.bg, curr_node)

            prev_end =  self.vectors[prev_node][side]
            prev_stem = (self.vectors[prev_node][(side+1)%2],self.vectors[prev_node][side])

            #print >>sys.stderr, "length:", length, "prev_stem:", prev_stem
            new_vec = add_close_segment(length, prev_stem)
            #print >>sys.stderr, "new_vec:", new_vec
            self.add_loop_segment(curr_node, prev_node, new_vec)

            self.vectors[curr_node] = new_vec
            self.sides[curr_node] = side

            #edge = list(bg.edges[curr_node])[0]

            for edge in self.bg.edges[curr_node]:
                if edge not in self.visited and edge != curr_node:
                    length = self.bg.get_stem_length(edge) * avg_stem_bp_length
                    #print >>sys.stderr, "length:", length
                    #print >>sys.stderr, "self.vectors:", self.vectors
                    segment = add_close_segment(length, new_vec)
                    self.vectors[edge] = segment
                    self.sides[edge] = side

                    self.add_stem_segment(edge, curr_node, segment)
                    self.visited.add(edge)

                    for edge1 in self.bg.edges[edge]:
                        if edge1 != edge:
                            self.to_visit.append((edge1, edge, side))

    def add_interior_loop(self, prev_node, curr_node, side):
        stem1 = self.vectors[prev_node]
        #print >>sys.stderr, "prev_node:", prev_node, "stem1:", stem1

        length = get_bulge_length(self.bg, curr_node, self.length_stats)
        bulge_dimensions = get_bulge_dimensions(self.bg, curr_node)

        #print >>sys.stderr, "bulge_dimensions:", bulge_dimensions
        dim0 = bulge_dimensions[0]
        dim1 = bulge_dimensions[1]

        found = False
        while not found:
            try:
                stats = choice(self.angle_stats[dim0][dim1])
                found = True
            except IndexError:
                dim1 -= 1

        angle1 = stats['angle_stem1_bulge']
        angle2 = stats['angle_stem2_bulge']
        angle3 = stats['angle_stem2_rotation']

        #print >>sys.stderr, "curr_node:", curr_node, "bulge_dimension[0]:", bulge_dimensions[0], "bulge_dimension[1]:", bulge_dimensions[1]
        #print >>sys.stderr, "angle1: %f, angle2: %f, angle3: %f" % (angle1, angle2, angle3)

        bulge = add_bulge(length, stem1, angle1)

        for edge in self.bg.edges[curr_node]:
            if edge != curr_node and edge not in self.visited:
                # the second stem
                l = self.bg.get_stem_length(edge) * avg_stem_bp_length
                stem2 = add_stem2(l, stem1, bulge, angle2, angle3)
                #print >>sys.stderr, "angle1: %f angle2: %f" % ( angle1, angle2 )
                self.vectors[edge] = stem2
                self.sides[edge] = side
                self.add_stem_segment(edge, curr_node, stem2)
                self.visited.add(edge)

                if self.bg.weights[curr_node] == 1:
                    self.add_single_bulge(curr_node, prev_node, bulge)

                if self.bg.weights[curr_node] == 2:
                    self.add_interior_bulge(curr_node, prev_node, bulge)

                self.vectors[curr_node] = bulge

                for edge1 in self.bg.edges[edge]:
                    # bulge
                    if edge1 != curr_node:
                        (sb1, se1) = self.bg.get_sides(edge, curr_node)
                        (sb2, se2) = self.bg.get_sides(edge, edge1)

                        if se1 == se2:
                            side = 0
                        else:
                            side = 1

                        self.to_visit.append((edge1, edge, side))
                        for edge2 in self.bg.edges[edge1]:
                            #stem
                            if edge2 in self.visited:
                                if list(stem2[0]) != list(self.vectors[edge2][0]):

                                    self.add_single_bulge(edge1, edge, (stem2[0], self.vectors[edge2][0]))

    def traverse_and_build(self):
        self.visited = set()
        self.vectors = dict()
        self.sides = dict()
        self.to_visit = []
        self.first_segment = True

        self.find_start_node()
        self.vectors['start'] = (array([0., 0., 0.]), array([0.,0.,1.]))

        counter = 0
        self.bg.coords = dict()

        while len(self.to_visit) > 0:
            (curr_node, prev_node, side) = self.to_visit.pop(0)

            while curr_node in self.visited:
                if len(self.to_visit) > 0:
                    (curr_node, prev_node, side) = self.to_visit.pop(0)
                else:
                    return
                    
            #print >>sys.stderr, "curr_node:", curr_node, "visited:", self.visited
            self.visited.add(curr_node)

            if curr_node[0] != 's':
                if len(self.bg.edges[curr_node]) == 1:
                    self.add_loop(prev_node, curr_node, side)

                if len(self.bg.edges[curr_node]) == 2:
                    self.add_interior_loop(prev_node, curr_node, side)

                elif self.bg.weights[curr_node] > 2:
                    self.add_something_else(prev_node, curr_node, side)


            #self.bg.output('building_%d.coords' % (counter))
            counter += 1
