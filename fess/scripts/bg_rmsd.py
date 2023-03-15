#!/usr/bin/python

from __future__ import absolute_import
from __future__ import print_function
import sys
from optparse import OptionParser

import tess.threedee.model.coarse_grain as ttmc
import borgy.graph.graph_pdb as cgg
import borgy.builder.rmsd as cbr

def main():
    usage = './bg_rmsd.py s1.coord s2.coord'
    parser = OptionParser()

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < 2:
        parser.print_help()
        sys.exit(1)

    bg1 = ttmc.CoarseGrainRNA(args[0])
    bg2 = ttmc.CoarseGrainRNA(args[1])

    r1 = cgg.bg_virtual_residues(bg1)
    r2 = cgg.bg_virtual_residues(bg2)

    print(cbr.centered_rmsd(r1, r2))

if __name__ == '__main__':
    main()

