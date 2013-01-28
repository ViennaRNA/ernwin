#!/usr/bin/python

import sys
import numpy as np
import warnings
import collections as co

import Bio.PDB as bp

import corgy.graph.bulge_graph as cgb
import corgy.graph.graph_pdb as cgg
import corgy.utilities.debug as cud

import corgy.visual.pymol as cvp

from optparse import OptionParser

'''
def add_pymol_coords(pp, coords):

    for i in range(2):
        for a in coords[i].keys():
            if a == 'P':
                pp.add_sphere(coords[i][a], 'orange', 0.1)
            elif a == 'C1':
                pp.add_sphere(coords[i][a], 'green', 0.1)
            elif a == 'O3*':
                pp.add_sphere(coords[i][a], 'red', 0.1)

def print_atom_positions(coords):
    for i in range(2):
        for a in coords[i].keys():
            print "%d %s %s" % (i, a,
                    " ".join(map(str, coords[i][a])))
'''


def print_average_atom_positions(coords, res='A', pp=None):
    averages = [co.defaultdict(lambda: co.defaultdict(list)), 
                co.defaultdict(lambda: co.defaultdict(list))]

    if pp == None:
        print "import collections as co"
        print "avg_stem_vres_atom_coords = [co.defaultdict(dict), co.defaultdict(dict)]\n"
    
    for i in range(2):
        for r in coords[i].keys():
            for c in coords[i][r]:
                for a in c.keys():
                    averages[i][r][a] += [c[a]]

    for i in range(2):
        for r in averages[i].keys():
            if r not in res:
                continue

            for a in averages[i][r].keys():
                if pp == None:
                    print "avg_stem_vres_atom_coords[%d]['%s']['%s'] = [%s]" % (i,r,a, 
                            ",".join(map(str, np.mean(averages[i][r][a], axis=0))))
                    #print i,r, a, np.mean(averages[i][r][a], axis=0)
                else:
                    if i == 0:
                        pp.add_sphere(np.mean(averages[i][r][a], axis=0), 'green', 0.1)
                    else:
                        pp.add_sphere(np.mean(averages[i][r][a], axis=0), 'red', 0.1)
        

def main():
    usage = './bounding_box_coords.py temp.comp temp.pdb'
    usage += "Print out the coordinates of two diagonal corners of the"
    usage += 'bounding box of each stem nucleotide in the pdb file.'
    parser = OptionParser()

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    parser.add_option('-e', '--edge', dest='edge', default=False, action='store_true', help='Include the edge nucleotides in the statistics.')
    parser.add_option('-p', '--pymol', dest='pymol', default=False, action='store_true', help='Output in pymol cgo format.')
    parser.add_option('-a', '--averages', dest='averages', default=False, action='store_true', help='Output the average coordinates of the bounding boxes')
    parser.add_option('-r', '--residue', dest='residue', default='A', help="The type of residue to calculate the averages for.", type='str')

    (options, args) = parser.parse_args()

    if len(args) < 2:
        parser.print_help()
        sys.exit(1)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        chain = list(bp.PDBParser().get_structure('temp', args[1]).get_chains())[0]

    bg = cgb.BulgeGraph(args[0])

    pp = cvp.PymolPrinter()
    pp.draw_axes = True

    all_coords = [co.defaultdict(list),
                  co.defaultdict(list)]

    for s in bg.stems():
        for i in range(bg.stem_length(s)):
            if i == 0 or i == bg.stem_length(s) - 1:
                if options.edge:
                    pass
                else:
                    continue

            (origin, basis, coords) = cgg.stem_vres_reference_atoms(bg, chain, s, i)

            # subtract one because the sequence is 0-based
            (p1, p2) = (bg.defines[s][0] + i - 1, bg.defines[s][3] - i - 1)
            (r1, r2) = (bg.seq[p1], bg.seq[p2])

            all_coords[0][r1] += [coords[0]]
            all_coords[1][r2] += [coords[1]]

            '''
            if options.pymol:
                if not options.averages:
                    add_pymol_coords(coords)
            else:
                print_atom_positions(coords)
            '''

    if options.averages:
        if options.pymol:
            print_average_atom_positions(all_coords, list(options.residue), pp)
            pp.output_pymol_file()
        else:
            print_average_atom_positions(all_coords, list(options.residue), None)
    else:
        for i, c in enumerate(all_coords):
            for k in c.keys():
                for coords in c[k]:
                    for atom in coords.keys():
                        print i, k, atom, " ".join(map(str, coords[atom]))


if __name__ == '__main__':
    main()

