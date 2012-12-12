#!/usr/bin/python

from optparse import OptionParser
import sys, pickle

import corgy.builder.energy as cbe

from corgy.graph.bulge_graph import BulgeGraph
from corgy.visual.pymol import PymolPrinter

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
    parser.add_option('-i', '--img_energy', dest='img_energy', default=False, action='store_true', help="Visualize the distance energy")
    parser.add_option('-s', '--stem_stem_energy', dest='stem_stem_energy', default=False, action='store_true', help="Visualize the distance energy")
    parser.add_option('-m', '--max_stem_distances', dest='max_stem_distances', default=0, help='Draw the vectors between the closest points on two different stems', type='float')
    parser.add_option('-p', '--pdb', dest='pdb_file', default=None, help='Include a pdb file for drawing bouding boxes.', type='string')
    parser.add_option('', '--hide-cg', dest='hide_cg', default=False, action='store_true', help='Hide the coarse grain model.')
    parser.add_option('', '--stem-atoms', dest='stem_atoms', default=False, action='store_true', help='Display the approximate locations of the atoms of the nucleotides that are parts of stems.')

    (options, args) = parser.parse_args()
    
    bgs = []
    for arg in args:
        bgs += [BulgeGraph(arg)]

    pymol_printer = PymolPrinter()
    pymol_printer.print_text = options.print_text
    pymol_printer.add_twists = options.add_twists
    pymol_printer.add_longrange = options.add_longrange
    pymol_printer.max_stem_distances = options.max_stem_distances

    if len(options.energy) > 0:
        for bg in bgs:
            bg.calc_bp_distances()
        pymol_printer.energy_function = pickle.load(open(options.energy, 'r'))

    if options.img_energy:
        pymol_printer.energy_function = cbe.ImgHelixOrientationEnergy()
    if options.stem_stem_energy:
        pymol_printer.energy_function = cbe.StemStemOrientationEnergy()
    if options.pdb_file:
        pymol_printer.pdb_file = options.pdb_file
    if options.hide_cg:
        pymol_printer.draw_segments = False

    for i in range(len(bgs)):
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

    if options.stem_atoms:
        for (s,i) in bg.virtual_residues():
            cud.pv('(s,i)')
            vra = cgg.virtual_residue_atoms(bg, s, i)
            (vpos1, vvec1) = cgg.virtual_res_3d_pos(bg, s, i)
            for i in range(2):
                for a in vra[i].values():
                    pymol_printer.add_sphere(a, 'purple', 0.3)

            for (s2, i2) in bg.virtual_residues():
                vra2 = cgg.virtual_residue_atoms(bg, s2, i2)
                (vpos2, vvec2) = cgg.virtual_res_3d_pos(bg, s2, i2)

                if s == s2:
                    continue

                for atoms1 in vra:
                    for atoms2 in vra2:
                        for a1 in atoms1.values():
                            for a2 in atoms2.values():
                                if cuv.magnitude(a1 - a2) < 2.0:
                                    mult = 7.
                                    cud.pv('(s,i,s2,i2)')
                                    cud.pv('cuv.magnitude((vpos1 + mult * vvec1) - (vpos2 + mult * vvec2))')
                                    pymol_printer.add_segment(a1, a2, 'yellow', 0.5)


    pymol_printer.output_pymol_file()

if __name__ == '__main__':
    main()
