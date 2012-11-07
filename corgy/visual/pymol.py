#!/usr/bin/python

import sys

import numpy as np

import corgy.graph.graph_pdb as cgg
import corgy.utilities.debug as cud
import corgy.utilities.vector as cuv

import Bio.PDB.Model as bpm
import Bio.PDB.Chain as bpc
import Bio.PDB.Structure as bps
import Bio.PDB as bp

class PymolPrinter:
    def __init__(self):
        self.new_segments = []
        self.segments = []
        self.spheres = []
        self.new_spheres = []
        self.override_color = None
        self.print_text = True
        self.energy_function = None
        self.add_twists = True
        self.add_longrange = False
        self.chain = None

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

    def add_sphere(self, p, color='green', width=0.2, text="", color_rgb = None):
        if self.override_color != None:
            color = self.override_color
        
        if color_rgb == None:
            color_rgb = self.get_color_vec(color)

        self.new_spheres += [(np.array(p), color_rgb, width, text)]


    def transform_spheres(self, translation, rotation):
        for (p, color, width, text) in self.new_spheres:
            p -= translation

            new_p = np.dot(rotation, p)
            self.spheres += [(p, color, width, text)]

        self.new_spheres = []

    def add_segment(self, p, n, color='green', width=0.2, text=""):
        if self.override_color != None:
            color = self.override_color

        #assert(not allclose(p, n))

        self.new_segments += [(np.array(p), np.array(n), color, width, text)]

    def transform_segments(self, translation, rotation):
        for (p, n, color, width, text) in self.new_segments:
            p -= translation
            n -= translation

            new_p = np.dot(rotation, p)
            new_n = np.dot(rotation, n)

            self.segments += [(new_p, new_n, color, width, text)]

        self.new_segments = []

    def pymol_spheres_string(self):
        self.spheres += self.new_spheres
        s = ''

        for (p, color, width, text) in self.new_spheres:
            color_vec = color
            s += "COLOR, %s," % (",  ".join([str(c) for c in color_vec])) + '\n'
            s += "SPHERE, %s, %f," % (", ".join([str(pi) for pi in p]), width) + '\n'

        return s
    
    def pymol_axis_string(self):
        w = 0.42 # cylinder width 
        l = 40.0 # cylinder length
        h = 3.0 # cone hight
        d = w * 2.618 # cone base diameter
        s = ""
         
        s += "CYLINDER, 0.0, 0.0, 0.0,   %f, 0.0, 0.0, %f, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0," % (l, w)
        s += "CYLINDER, 0.0, 0.0, 0.0, 0.0,   %f, 0.0, %f, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0," % (l, w)
        s += "CYLINDER, 0.0, 0.0, 0.0, 0.0, 0.0,   %f, %f, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0," % (l, w)
        s += "CONE,   %f, 0.0, 0.0, %f, 0.0, 0.0, %f, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0," % (l, h+l, d)
        s += "CONE, 0.0, %f, 0.0, 0.0, %f, 0.0, %f, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0," % (l, h+l, d)
        s += "CONE, 0.0, 0.0, %f, 0.0, 0.0, %f, %f, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0, 1.0," % (l, h+l, d)

        return s

    def pymol_segments_string(self):
        color = 'green'
        width = 0.2
        s = ''

        self.segments += self.new_segments

        for seg  in self.segments:
            (p,n,color,width, text) = seg
            color_vec = [str(c) for c in self.get_color_vec(color)]
            s += " CYLINDER, %f, %f, %f, %f, %f, %f, %f, %s, %s," % (p[0], p[1], p[2], n[0], n[1], n[2], width, ", ".join(color_vec), ", ".join(color_vec)) + '\n'

        return s

    def pymol_text_string(self):
        counter = 0
        s = ''

        for (p, n, color, width, text) in self.segments:
            if len(text) == 0:
                continue

            s +=  "cgox_%d = []" % (counter) + '\n'

            comp1 = cuv.normalize(n - p)

            ncl = cuv.get_non_colinear_unit_vector(comp1)

            comp2 = cuv.normalize(np.cross(ncl, comp1))
            comp3 = cuv.normalize(np.cross(ncl, comp2))

            #pos = (p + n) / 2.0 + 3 * comp2
            pos = p + (n - p) / 4.0 + 3 * comp2
            font = 2
            axes = [list(comp1 * 2), list(comp2 * 2), list(comp3 * 2)]

            text = "%s: %.1f" % (text, cuv.magnitude(n-p))

            s += "cyl_text(cgox_%d, plain, %s, \"%s\", 0.20, axes=%s)" % (counter, str(list(pos)), text, str(axes)) + '\n'
            counter += 1

        s +=  "cmd.set(\"cgo_line_radius\",0.03)" + '\n'
        for i in range(counter):
            s += "cmd.load_cgo(cgox_%d, \'cgox%d\')" % (i, i) + '\n'
        s += "cmd.zoom(\"all\", 2.0)" + '\n'

        return s


    def pymol_string(self):
        '''
        Output the contents of this structure into a file that can be passed
        in as a pymol script.
        '''

        s = self.pymol_intro_string()
        s += self.pymol_segments_string()
        s += self.pymol_spheres_string()
        s += self.pymol_axis_string()
        s += self.pymol_outro_string()

        if self.print_text:
            s += self.pymol_text_string()

        return s

    def dump_pdb(self, filename):
        '''
        If the BulgeGraph has a chain created for it, dump that as well.

        @param filename: The filename of the pdb file to which the chain coordinates will be written.
        '''
        if self.chain == None:
            return

        self.chain.child_list.sort()
        m = bpm.Model(' ')
        s = bps.Structure(' ')

        m.add(self.chain)
        s.add(m)

        io = bp.PDBIO()
        io.set_structure(s)
        io.save(filename)


    def dump_pymol_file(self, filename):
        '''
        Output the structure to file.

        @param filename: The location of the output file.
        '''
        # Output the script for showing the coarse-grained elements
        f = open(filename + ".pym", 'w')
        f.write(self.pymol_string())
        f.close()

        # Output the pdb structure
        self.dump_pdb(filename + ".pdb")

        # Output the script file for loading the pdb and coarse grained structure
        f = open(filename + ".pml", 'w')
        f.write("run %s" % (filename + ".pym"))
        f.close()

    def output_pymol_file(self):
        print self.pymol_string()

    def reset(self):
        self.segments = []
        self.new_segments = []

    def pymol_intro_string(self):
        s  = "from pymol.cgo import *" + '\n'
        s += "from pymol import cmd" + '\n'
        s += "from pymol.vfont import plain" + '\n'
        s += "obj = [" + '\n'
        return s

    def pymol_outro_string(self):
        s =  "]" + '\n'
        s += "cmd.load_cgo(obj, 'ss')" + '\n'
        return s

    def add_stem_like(self, bg, key, color = 'green', width=2.4):
        (p, n) = bg.coords[key]
        self.add_segment(p, n, color, width, key)

        if self.add_twists:
            twist1o = bg.get_twists(key)[0]
            twist2o = bg.get_twists(key)[1]

            twist_rot_mat_l = cuv.rotation_matrix(n - p, -(1.45 / 2.))
            twist_rot_mat_r = cuv.rotation_matrix(n - p, (1.45 / 2.))

            twist1 = np.dot(twist_rot_mat_l, twist1o)
            twist2 = np.dot(twist_rot_mat_l, twist2o)

            twist3 = np.dot(twist_rot_mat_r, twist1o)
            twist4 = np.dot(twist_rot_mat_r, twist2o)


            mult = 7.
            width = .3

            self.add_segment(p, p + mult * twist1, "white", width, '')
            self.add_segment(n, n + mult * twist2, "white", width, '')

            self.add_segment(p, p + mult * twist3, "red", width, '')
            self.add_segment(n, n + mult * twist4, "red", width, '')

        stem_len = bg.defines[key][1] - bg.defines[key][0] + 1
        for i in range(stem_len):
            (pos, vec) = cgg.virtual_res_3d_pos(bg, key, i)
            self.add_segment(pos, pos + mult * vec, "blue", width, '')

        '''
        self.add_sphere(p + mult * twist1, "white", width, key)
        self.add_sphere(n + mult * twist2, "white", width, key)
        '''

    def coordinates_to_pymol(self, bg):
        for key in bg.coords.keys():
            (p, n) = bg.coords[key]
        
            if key[0] == 's':
                self.add_stem_like(bg, key)

            else:
                if len(bg.edges[key]) == 1:
                    self.add_segment(p, n, "blue", 1.0, key)
                if len(bg.edges[key]) == 2:
                    if bg.weights[key] == 1:
                        self.add_segment(p, n, "red", 1.0, key + " " + str(bg.defines[key][1] - bg.defines[key][0]) + "")
                    else:
                        self.add_stem_like(bg, key, "yellow", 1.0)
                        #self.add_segment(p, n, "yellow", 1.0, key)

        if self.add_longrange:
            for key1 in bg.longrange.keys():
                for key2 in bg.longrange[key1]:
                    try:
                        point1 = bg.get_point(key1)
                        point2 = bg.get_point(key2)

                        self.add_segment(point1, point2, "purple", 0.3, key1 + " " + key2)
                        self.add_segment(point1, point2, "purple", 0.3, key1 + " " + key2)
                    except:
                        continue

        print >>sys.stderr, "energy_function:", self.energy_function
        # print the contributions of the energy function, if one is specified
        if self.energy_function != None:
            print >>sys.stderr, "key"
            sum_energy = 0.

            int_energies = list(self.energy_function.iterate_over_interaction_energies(bg, background=False))
            min_energy = min(int_energies, key=lambda x: x[1])
            print >>sys.stderr, "min_energy:", min_energy

            for (interaction, energy) in int_energies:
                (p, n) = (bg.get_point(interaction[0]), bg.get_point(interaction[1]))
                scaled_energy = min_energy[1] - energy

                print >>sys.stderr, "Adding segment:", energy, scaled_energy, np.exp(scaled_energy) 
                self.add_segment(p, n, 'purple', 3 * np.exp(scaled_energy) )

                sum_energy += energy

            cud.pv("sum_energy")

    def chain_to_pymol(self, chain):
        '''
        Add a Bio.PDB.Chain to the structure, so that it can later be printed.
        '''
        self.chain = chain

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

