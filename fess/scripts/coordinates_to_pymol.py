#!/usr/bin/python

from optparse import OptionParser
import sys, pickle

import corgy.builder.energy as cbe
import numpy as np

from corgy.graph.bulge_graph import BulgeGraph
from corgy.visual.pymol import PymolPrinter

import corgy.builder.models as cbm
import corgy.graph.graph_pdb as cgg
import corgy.utilities.debug as cud
import corgy.utilities.vector as cuv
import corgy.utilities.average_stem_vres_atom_positions as cua

def main():
    if len(sys.argv) < 2:
        print "Usage: ./coordinates_to_pymol temp.coordinates"
        sys.exit(1)

    parser = OptionParser()

    parser.add_option('-c', '--centers', dest='centers', default=0, help='Display the centers of each segment.', type='int') 
    parser.add_option('-f', '--flexibility', dest='flex', default=None, help='Location of the flexibility statistics file.', type='string') 
    parser.add_option('-x', '--text', dest='print_text', default=False, action='store_true', help='Print the names of the segments in the pymol output')
    parser.add_option('-t', '--twists', dest='add_twists', default=True, action='store_false', help='Hide the twist indicators')
    parser.add_option('-l', '--longrange', dest='add_longrange', default=False, action='store_true', help='Display the longrange interactions')
    parser.add_option('-e', '--energy', dest='energy', default='', help='Location of an energy function to visualize', type='string')
    parser.add_option('-i', '--img-energy', dest='img_energy', default=False, action='store_true', help="Visualize the distance energy")
    parser.add_option('-s', '--stem-stem-energy', dest='stem_stem_energy', default=False, action='store_true', help="Visualize the distance energy")
    parser.add_option('', '--cylinder-intersect-energy', dest='cylinder_intersect_energy', default=False, action='store_true', help="Visualize the distance energy")
    parser.add_option('-m', '--max-stem-distances', dest='max_stem_distances', default=0, help='Draw the vectors between the closest points on two different stems', type='float')
    parser.add_option('-p', '--pdb', dest='pdb_file', default=None, help='Include a pdb file for drawing bouding boxes.', type='string')
    parser.add_option('', '--hide-cg', dest='hide_cg', default=False, action='store_true', help='Hide the coarse grain model.')
    parser.add_option('', '--stem-atoms', dest='stem_atoms', default=False, action='store_true', help='Display the approximate locations of the atoms of the nucleotides that are parts of stems.')
    parser.add_option('', '--stem-atom-distances', dest='stem_atom_distances', default=False, action='store_true', help='Check the distances between the stem atoms.')
    parser.add_option('', '--movie', dest='movie', default=False, action='store_true')
    parser.add_option('', '--align-to-longest-stem', dest='align_to_longest_stem', default=False, action='store_true', help='Align all of the models to the longest stem of the first one')
    parser.add_option('', '--include-pdb', dest='include_pdb', default=False, action='store_true', help='Include the pdb files in a movie. The pdb files should be named just like the graph files except with a .pdb appended.')

    (options, args) = parser.parse_args()
    
    bgs = []
    for arg in args:
        bgs += [BulgeGraph(arg)]

    pymol_printer = PymolPrinter()
    pymol_printer.print_text = options.print_text
    pymol_printer.add_twists = options.add_twists
    pymol_printer.add_longrange = options.add_longrange
    pymol_printer.max_stem_distances = options.max_stem_distances
    pymol_printer.movie = options.movie

    if len(options.energy) > 0:
        for bg in bgs:
            bg.calc_bp_distances()
        pymol_printer.energy_function = pickle.load(open(options.energy, 'r'))

    if options.img_energy:
        pymol_printer.energy_function = cbe.ImgHelixOrientationEnergy()
    if options.stem_stem_energy:
        for bg in bgs:
            for s in bg.stems():
                cgg.add_virtual_residues(bg, s)
        pymol_printer.energy_function = cbe.StemStemOrientationEnergy()
    if options.cylinder_intersect_energy:
        for bg in bgs:
            for s in bg.stems():
                cgg.add_virtual_residues(bg, s)
        pymol_printer.energy_function = cbe.CylinderIntersectionEnergy()
    if options.pdb_file:
        pymol_printer.pdb_file = options.pdb_file
    if options.hide_cg:
        pymol_printer.draw_segments = False

    for i in range(len(bgs)):
        bg = bgs[i]
        pymol_printer.reset()

        if not options.movie:
            if len(bgs) > 1 and i == 0:
                pymol_printer.override_color = 'green'
            elif len(bgs) > 1 and i > 0:
                pymol_printer.override_color = 'red'

        #pymol_printer.override_color = 'red'
        bg = bgs[i]

        if options.centers == 1:
            pymol_printer.centers_to_pymol(bg)
        elif options.flex != None:
            pymol_printer.flex_to_pymol(bg, options.flex)
        
        pymol_printer.coordinates_to_pymol(bg)

        if options.align_to_longest_stem:
            sm = cbm.SpatialModel(bg)
            longest_stem = bgs[0].get_longest_stem()

            sm.stems[longest_stem] = cbm.StemModel(name=longest_stem)
            sm.stems[longest_stem].mids = bg.coords[longest_stem]
            sm.stems[longest_stem].twists = bg.twists[longest_stem]

            translation_matrix = sm.get_transform(longest_stem)
            rotation_matrix = sm.get_rotation(longest_stem)
        else:
            translation_matrix = np.zeros(3)
            rotation_matrix = np.identity(3)

        #cud.pv('translation_matrix')
        #cud.pv('rotation_matrix')
        pymol_printer.transform_segments(translation_matrix, rotation_matrix)

        if options.stem_atoms:
            for s in bg.stems():
                cgg.add_virtual_residues(bg, s)

            for (s,i) in bg.virtual_residues():
                vra = []
                vra += [cgg.virtual_residue_atoms(bg, s, i,0)]
                vra += [cgg.virtual_residue_atoms(bg, s, i,1)]


                #(vpos1, vvec1, vvec1_l, vvec1_r) = cgg.virtual_res_3d_pos(bg, s, i)
                (vpos1, vvec1, vvec1_l, vvec1_r) = bg.v3dposs[s][i]

                '''
                print >> sys.stderr, "************************"
                if s == 's1' and i == 1:
                    print >> sys.stderr, "----------------------"
                    cud.pv('vpos1')
                    cud.pv('vvec1_l')
                    cud.pv('vpos1 + mult * vvec1_l')
                '''

                for x in range(2):
                    for a in vra[x].values():
                        pymol_printer.add_sphere(a, 'purple', 0.3)

                if not options.stem_atom_distances:
                    continue

                #cud.pv('(s,i, s2, i2)')
                for (s2, i2) in bg.virtual_residues():
                    vra2 = []
                    vra2 += [cgg.virtual_residue_atoms(bg, s2, i2,0)]
                    vra2 += [cgg.virtual_residue_atoms(bg, s2, i2,1)]
                    #(vpos2, vvec2, vvec2_l, vvec2_r) = cgg.virtual_res_3d_pos(bg, s2, i2)
                    (vpos2, vvec2, vvec2_l, vvec2_r) = bg.v3dposs[s2][i2]

                    if s == s2:
                        continue

                    #for atoms1 in vra:
                        #for atoms2 in vra2:
                    for k in range(2):
                        atoms1 = vra[k]
                        for m in range(2):
                            atoms2 = vra2[m]
                            for a1 in atoms1.values():
                                for a2 in atoms2.values():
                                    if cuv.magnitude(a1 - a2) < 2.0:
                                        mult = 8.
                                        cud.pv('(s,i,k,s2,i2,m)')
                                        '''
                                        cud.pv('(k,m)')

                                        cud.pv('vpos1')
                                        cud.pv('vvec1_l')
                                        cud.pv('vpos1 + mult * vvec1_l')
                                        cud.pv('vpos1 + mult * vvec1_r')
                                        cud.pv('vpos2 + mult * vvec2_l')
                                        cud.pv('vpos2 + mult * vvec2_r')
                                        '''
                                        cud.pv('cuv.magnitude((vpos1 + mult * vvec1_l) - (vpos2 + mult * vvec2_l))')
                                        cud.pv('cuv.magnitude((vpos1 + mult * vvec1_l) - (vpos2 + mult * vvec2_r))')
                                        cud.pv('cuv.magnitude((vpos1 + mult * vvec1_r) - (vpos2 + mult * vvec2_l))')
                                        cud.pv('cuv.magnitude((vpos1 + mult * vvec1_r) - (vpos2 + mult * vvec2_r))')
                                        cud.pv('cuv.magnitude(a1 - a2)')
                                        pymol_printer.add_segment(a1, a2, 'yellow', 0.5)


        pymol_printer.output_pymol_file()

        if options.include_pdb:
            print 'cmd.load("%s", "sxs%s", state=%d)' % (args[i] + ".pdb", pymol_printer.prev_obj_name, pymol_printer.state-1)
            print 'cmd.color("grey50", "all")'

    

if __name__ == '__main__':
    main()
