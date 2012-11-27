#!/usr/bin/python

import sys
from optparse import OptionParser

import corgy.graph.graph_pdb as cgg
import corgy.visual.pymol as cvp

def print_base_normals(pdb_file):
    '''
    Calculate the normal of each base in the pdb file
    and output a line through its center normal to its
    place as defined by the C2-N3 and C2-C5 vectors.
    '''
    base_normals = cgg.base_normals(pdb_file)
    pp = cvp.PymolPrinter()

    for (center, norm) in base_normals:
        norm_len = 14.

        p = center + norm_len * norm
        n = center - norm_len * norm

        pp.add_segment(p, n, "white", 0.5, '')

    pp.output_pymol_file()

def main():
    usage = './base_normals_to_pymol.py pdb_file\n'
    usage += 'Display cylinders where the base normals should be.'

    parser = OptionParser()

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    pdb_file = args[0]
    print >>sys.stderr, "pdb_file", pdb_file
    print_base_normals(pdb_file)

if __name__ == '__main__':
    main()

