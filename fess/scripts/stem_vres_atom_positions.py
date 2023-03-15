#!/usr/bin/python

from __future__ import absolute_import
from __future__ import print_function
import sys
import numpy as np
import warnings
import collections as co
import itertools

import Bio.PDB as bp

import forgi.threedee.model.coarse_grain as ttmc
import forgi.threedee.utilities.graph_pdb as cgg
import forgi.utilities.debug as cud

import forgi.threedee.visual.pymol as cvp

from optparse import OptionParser
import matplotlib.pyplot as plt

import logging
from six.moves import map
from six.moves import range
logging.basicConfig(level=logging.WARNING)
log = logging.getLogger(__name__)
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
    
    if pp is None:
        print("import collections as co")
        print("avg_stem_vres_atom_coords = [co.defaultdict(dict), co.defaultdict(dict)]\n")
    
    for i in range(2):
        for k in coords[i].keys():
            for c in coords[i][k]:
                for a in c.keys():
                    averages[i][k][a] += [c[a]]


    import matplotlib.colors
    fig, ax = plt.subplots(3)
    colors = list(matplotlib.colors.cnames.keys())
    markers = itertools.cycle(["+", "o", ".", "<", ">", "^", "v", "s", "x", "d", "D", "*"])
    for i in range(2):
        for r in averages[i].keys():
            if r not in res:
                continue

            for j, a in enumerate(averages[i][r].keys()):
                if pp is None:
                    print('avg_stem_vres_atom_coords[%d]["%s"]["%s"] = [%s]' % (i,r,a, 
                            ",".join(map(str, np.mean(averages[i][r][a], axis=0)))))
                    #print i,r, a, np.mean(averages[i][r][a], axis=0)
                else:
                    if i == 0:
                        pp.add_sphere(np.mean(averages[i][r][a], axis=0), 'green', 0.1)
                    else:
                        pp.add_sphere(np.mean(averages[i][r][a], axis=0), 'red', 0.1)
                    pp.add_text(np.mean(averages[i][r][a], axis=0), a, color="yellow" )
                if i==0:
                    ax[0].scatter(np.array(averages[i][r][a])[:,0], np.array(averages[i][r][a])[:,1], label="{}:{}".format(i, a), color=colors[j], marker = next(markers))
                    ax[1].scatter(np.array(averages[i][r][a])[:,1], np.array(averages[i][r][a])[:,2], label="{}:{}".format(i, a), color=colors[j], marker = next(markers))
                    ax[2].scatter(np.array(averages[i][r][a])[:,0], np.array(averages[i][r][a])[:,2], label="{}:{}".format(i, a), color=colors[j], marker = next(markers))

    ax[0].set_xlabel("x")
    ax[0].set_ylabel("y")
    ax[1].set_xlabel("y")
    ax[1].set_ylabel("z")
    ax[2].set_xlabel("x")
    ax[2].set_ylabel("z")

    ax[2].legend()
    plt.show()

def main():
    usage = './bounding_box_coords.py temp.pdb [temp2.pdb ...]'
    usage += "Print out the positions of the atoms in coordinates "
    usage += "respective to the virtual residue coordinate system."
    parser = OptionParser()

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    parser.add_option('-e', '--edge', dest='edge', default=False, action='store_true', help='Include the edge nucleotides in the statistics.')
    parser.add_option('-p', '--pymol', dest='pymol', default=False, action='store_true', help='Output in pymol cgo format.')
    parser.add_option('-a', '--averages', dest='averages', default=False, action='store_true', help='Output the average coordinates')
    parser.add_option('-r', '--residue', dest='residue', default='AUGC', help="The type of residue to calculate the averages for.", type='str')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    all_coords = [co.defaultdict(list),
                      co.defaultdict(list)]

    pp = cvp.PymolPrinter()
    pp.draw_axes = True

    for pdbfile in args:
        try:
            bg = ttmc.from_pdb(pdbfile)
        except Exception as e:
            log.exception(e)
            continue 
        chain = bg.chain

        for s in bg.stem_iterator():
            for i in range(bg.stem_length(s)):
                if i == 0 or i == bg.stem_length(s) - 1:
                    if not options.edge:
                        continue

                (origin, basis, coords) = cgg.stem_vres_reference_atoms(bg, chain, s, i)

                # subtract one because the sequence is 0-based
                (p1, p2) = (bg.defines[s][0] + i - 1, bg.defines[s][3] - i - 1)
                (r1, r2) = (bg.seq[p1], bg.seq[p2])

                all_coords[0][r1].append(coords[0])
                all_coords[1][r2].append(coords[1])

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
                        print(i, k, atom, " ".join(map(str, coords[atom])))


if __name__ == '__main__':
    main()

