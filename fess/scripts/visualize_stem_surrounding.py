#!/usr/bin/python

import sys, os
from optparse import OptionParser

import numpy as np

import corgy.graph.bulge_graph as cgb
import corgy.graph.graph_pdb as cgg
import corgy.builder.config as cbc
import corgy.visual.pymol as cvp
import corgy.utilities.vector as cuv

def main():
    usage = """
Visualize the data about the location of surrounding virutal
residues.

usage: %prog [options] data_file
"""
    parser = OptionParser(usage = usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    parser.add_option('-w', '--with_stem', dest='with_stem', default=False, action='store_true', help='Draw the surrounding residue distribution with respect to a stem')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        print >>sys.stderr, "Lack of data_file argument."
        parser.print_help()
        sys.exit(1)

    if len(args) > 3:
        print >>sys.stderr, "Too many files. At most 3 can be displayed at once."
        sys.exit(1)


    pp = cvp.PymolPrinter()
    colors = ['red', 'green', 'blue']

    for i in range(len(args)):
        f = open(args[i], 'r')
        for line in f.readlines():
            parts = line.split(' ')
            res_type = int(parts[0])
            x = float(parts[2])
            y = float(parts[3])
            z = float(parts[4])
            
            pp.add_sphere([x,y,z], colors[i], 0.3, '')

    pp.output_pymol_file()

if __name__ == '__main__':
    main()

