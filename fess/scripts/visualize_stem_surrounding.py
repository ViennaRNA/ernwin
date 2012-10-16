#!/usr/bin/python

import sys, os
from optparse import OptionParser

import numpy as np
import scipy.stats as stats

import corgy.graph.bulge_graph as cgb
import corgy.graph.graph_pdb as cgg
import corgy.builder.config as cbc
import corgy.visual.pymol as cvp
import corgy.utilities.vector as cuv
import corgy.utilities.colormap as cuc
import corgy.exp.kde as cek

import pandas as pa

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

    kernels = []
    point_sets = []
    val_sets = []
    kernels = []

    for i in range(len(args)):
        stats = pa.read_csv(args[i],header=None, sep=' ')
        points = stats[['X.3', 'X.4', 'X.5']].as_matrix()

        point_sets += [points]

        kernel = cek.gaussian_kde(points.T)
        kernels += [kernel]

        vals = kernel(points.T)
        val_sets += [vals]


    if len(val_sets) == 2:
        vals = np.log(kernels[0](point_sets[0].T)) - np.log(kernels[1](point_sets[0].T))
    else:
        vals = np.log(kernels[0](point_sets[0].T))

    mi = min(vals)
    mx = max(vals)

    for k in range(1):
        for j in range(len(point_sets[0])):
            #print "p:", p, "kernel(p):", kernel(p)
            val = vals[j]
            pp.add_sphere(point_sets[k][j], colors[i], 0.3, '', color_rgb=cuc.floatRgb(val, mi, mx))

    pp.output_pymol_file()

if __name__ == '__main__':
    main()

