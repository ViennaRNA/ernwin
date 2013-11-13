#!/usr/bin/python

import sys
import numpy as np

import Bio.PDB as bp

import tess.threedee.model.coarse_grain as ttmc
import borgy.graph.graph_pdb as cgg

import borgy.visual.pymol as cvp

from optparse import OptionParser

def main():
    usage = './bounding_box_coords.py temp.comp temp.pdb'
    usage += "Print out the coordinates of two diagonal corners of the"
    usage += 'bounding box of each stem nucleotide in the pdb file.'
    parser = OptionParser()

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    parser.add_option('-e', '--edge', dest='edge', default=False, action='store_true', help='Include the edge nucleotides in the statistics.')
    parser.add_option('-p', '--pymol', dest='pymol', default=False, action='store_true', help='Output in pymol cgo format.')
    parser.add_option('-a', '--averages', dest='averages', default=False, action='store_true', help='Output the average coordinates of the bounding boxes')

    (options, args) = parser.parse_args()

    if len(args) < 2:
        parser.print_help()
        sys.exit(1)

    chain = list(bp.PDBParser().get_structure('temp', args[1]).get_chains())[0]
    bg = ttmc.CoarseGrainRNA(args[0])

    pp = cvp.PymolPrinter()
    pp.draw_axes = True
    all_corners = []

    for s in bg.stems():
        for i in range(bg.stem_length(s)):
            if i == 0 or i == bg.stem_length(s) - 1:
                if options.edge:
                    pass
                else:
                    continue

            (origin, basis, corners) = cgg.bounding_boxes(bg, chain, s, i)
            all_corners += [corners]

            if options.pymol:
                pp.add_sphere(corners[0][0], 'yellow', 0.1, '', [238/255., 221/255., 130/255.])
                pp.add_sphere(corners[0][1], 'yellow', 0.1, '', [184/255.,134/255.,11/255.])

                pp.add_sphere(corners[1][0], 'purple', 0.1, '', [238/255., 130/255., 238/255.])
                pp.add_sphere(corners[1][1], 'purple', 0.1, '', [208/255., 32/255., 144/255.])
            else:
                print '0','0', " ".join(map(str, corners[0][0]))
                print '0','1', " ".join(map(str, corners[0][1]))
                print '1','0', " ".join(map(str, corners[1][0]))
                print '1','1', " ".join(map(str, corners[1][1]))

    if options.averages:
        all_corners = np.array(all_corners)
        print '--averages--'
        print '0', '0', " ".join(map(str, np.mean(np.array([c[0][0] for c in all_corners]), axis=0)))
        print '0', '1', " ".join(map(str, np.mean(np.array([c[0][1] for c in all_corners]), axis=0)))
        print '1', '0', " ".join(map(str, np.mean(np.array([c[1][0] for c in all_corners]), axis=0)))
        print '1', '1', " ".join(map(str, np.mean(np.array([c[1][1] for c in all_corners]), axis=0)))

    if options.pymol:
        pp.output_pymol_file()

if __name__ == '__main__':
    main()

