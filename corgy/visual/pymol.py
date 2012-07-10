#!/usr/bin/python

from numpy import array, cross, dot, allclose, exp
from corgy.utilities.vector import normalize, get_non_colinear_unit_vector, magnitude
from sys import stderr

class PymolPrinter:
    def __init__(self):
        self.new_segments = []
        self.segments = []
        self.spheres = []
        self.new_spheres = []
        self.override_color = None
        self.print_text = True
        self.energy_function = None

    def get_color_vec(self, color):
        if color == 'green':
            return [0.0, 1.0, 0.0]
        elif color == 'blue':
            return [0.0, 0.0, 1.0]
        elif color == 'red':
            return [1.0, 0.0, 0.0]
        elif color == 'yellow':
            return [1.0, 1.0, 0.0]
        elif color == 'purple':
            return [1.0, 0.0, 1.0]
        elif color == 'white':
            return [1.0, 1.0, 1.0]
        else:
            return [0.0, 0.0, 0.0]

    def add_sphere(self, p, color='green', width=0.2, text=""):
        if self.override_color != None:
            color = self.override_color

        self.new_spheres += [(array(p), color, width, text)]


    def transform_spheres(self, translation, rotation):
        for (p, color, width, text) in self.new_spheres:
            p -= translation

            new_p = dot(rotation, p)
            self.spheres += [(p, color, width, text)]

        self.new_spheres = []

    def add_segment(self, p, n, color='green', width=0.2, text=""):
        if self.override_color != None:
            color = self.override_color

        #assert(not allclose(p, n))

        self.new_segments += [(array(p), array(n), color, width, text)]

    def transform_segments(self, translation, rotation):
        for (p, n, color, width, text) in self.new_segments:
            p -= translation
            n -= translation

            new_p = dot(rotation, p)
            new_n = dot(rotation, n)

            self.segments += [(new_p, new_n, color, width, text)]

        self.new_segments = []

    def print_pymol_spheres(self):
        self.spheres += self.new_spheres

        for (p, color, width, text) in self.new_spheres:
            color_vec = self.get_color_vec(color)
            print "COLOR, %s," % (",  ".join([str(c) for c in color_vec]))
            print "SPHERE, %s, %f," % (", ".join([str(pi) for pi in p]), width)

    def print_pymol_segments(self):
        color = 'green'
        width = 0.2

        self.segments += self.new_segments

        for seg  in self.segments:
            (p,n,color,width, text) = seg
            color_vec = [str(c) for c in self.get_color_vec(color)]
            print " CYLINDER, %f, %f, %f, %f, %f, %f, %f, %s, %s," % (p[0], p[1], p[2], n[0], n[1], n[2], width, ", ".join(color_vec), ", ".join(color_vec))


    def print_pymol_text(self):
        counter = 0

        for (p, n, color, width, text) in self.segments:
            if len(text) == 0:
                continue

            print "cgox_%d = []" % (counter)

            comp1 = normalize(n - p)

            ncl = get_non_colinear_unit_vector(comp1)

            comp2 = normalize(cross(ncl, comp1))
            comp3 = normalize(cross(ncl, comp2))

            pos = (p + n) / 2.0 + 3 * comp2
            font = 2
            axes = [list(comp1 * 2), list(comp2 * 2), list(comp3 * 2)]

            text = "%s: %.1f" % (text, magnitude(n-p))

            print "cyl_text(cgox_%d, plain, %s, \"%s\", 0.20, axes=%s)" % (counter, str(list(pos)), text, str(axes))
            counter += 1

        print "cmd.set(\"cgo_line_radius\",0.03)"
        for i in range(counter):
            print "cmd.load_cgo(cgox_%d, \'cgox%d\')" % (i, i)
        print "cmd.zoom(\"all\", 2.0)"

    def output_pymol_file(self):
        self.print_pymol_intro()
        self.print_pymol_segments()
        self.print_pymol_spheres()
        self.print_pymol_outro()

        if self.print_text:
            self.print_pymol_text()

    def reset(self):
        self.segments = []
        self.new_segments = []

    def print_pymol_intro(self):
        print "from pymol.cgo import *"
        print "from pymol import cmd"
        print "from pymol.vfont import plain"
        print "obj = ["

    def print_pymol_outro(self):
        print "]"
        print "cmd.load_cgo(obj, 'ss')"

    def coordinates_to_pymol(self, bg):
        for key in bg.coords.keys():
            (p, n) = bg.coords[key]
        
            if key[0] == 's':
                self.add_segment(p, n, 'green', 2.4, key)

                twist1 = bg.twists[key][0]
                twist2 = bg.twists[key][1]


                mult = 5.
                width = .3

                self.add_segment(p, p + mult * twist1, "white", width, '')
                self.add_segment(n, n + mult * twist2, "white", width, '')

                '''
                self.add_sphere(p + mult * twist1, "white", width, key)
                self.add_sphere(n + mult * twist2, "white", width, key)
                '''

            else:
                if len(bg.edges[key]) == 1:
                    self.add_segment(p, n, "blue", 1.0, key)
                if len(bg.edges[key]) == 2:
                    if bg.weights[key] == 1:
                        self.add_segment(p, n, "red", 1.0, key)
                    else:
                        self.add_segment(p, n, "yellow", 1.0, key)

        for key1 in bg.longrange.keys():
            for key2 in bg.longrange[key1]:
                try:
                    point1 = bg.get_point(key1)
                    point2 = bg.get_point(key2)

                    #self.add_segment(point1, point2, "purple", 0.3, key1 + " " + key2)
                    #self.add_segment(point1, point2, "purple", 0.3, key1 + " " + key2)
                except:
                    continue

        # print the contributions of the energy function, if one is specified
        if self.energy_function != None:
            for (interaction, energy) in self.energy_function.iterate_over_interactions(bg, background=False):
                print >>stderr, interaction, energy
                (p, n) = (bg.get_point(interaction[0]), bg.get_point(interaction[1]))
                self.add_segment(p, n, 'purple', exp(energy) * 10)

    def load_flex_stats(self, flex_file):
        f = open(flex_file, 'r')

        d = DefaultDict(DefaultDict(0.))

        for line in f.readlines():
            parts = line.strip().split(' ')

            d[int(parts[0])][int(parts[1])] = float(parts[2])

        return d

    def flex_to_pymol(self, bg, flex_file):
        flex_stats = self.load_flex_stats(flex_file)

        for key in bg.defines.keys():
            if key[0] != 's':
                if key in bg.coords.keys():
                    coords = bg.coords[key]
                    p = (coords[1] + coords[0]) / 2.

                    bd = bg.defines[key]
                    
                    if len(bd) == 2:
                        #out_str += "0 %d" % (abs(bd[1] - bd[0]) + 1)
                        dims = (0, abs(bd[1] - bd[0]) + 1)
                    else:
                        dims = (abs(bd[1] - bd[0]) + 1, abs(bd[2] - bd[3]) + 1)
                        #out_str += "%d %d" % ( min(dims), max(dims))

                    flex = flex_stats[min(dims)][max(dims)] * 10.

                    if len(bg.edges[key]) == 2:
                        if bg.weights[key] == 1:
                            self.add_sphere(p, "red", flex, key)
                        else:
                            self.add_sphere(p, "yellow", flex, key)


    def centers_to_pymol(self, bg):
        for key in bg.defines.keys():
            if key in bg.coords.keys():
                coords = bg.coords[key]
                p = (coords[1] + coords[0]) / 2.

                if key[0] == 's':
                    self.add_sphere(p, 'green', 3, key)
                else:
                    if len(bg.edges[key]) == 1:
                        self.add_sphere(p, 'blue', 1.5, key)
                if len(bg.edges[key]) == 2:
                    if bg.weights[key] == 1:
                        self.add_sphere(p, "red", 1.5, key)
                    else:
                        self.add_sphere(p, "yellow", 1.5, key)


def print_angle_stats():
    angles = []

    for i in range(0, len(segments)):
        s1 = segments[i-1][1] - segments[i-1][0]
        s2 = segments[i][1] - segments[i][0]

        angles += [vec_angle(s1, s2)]

