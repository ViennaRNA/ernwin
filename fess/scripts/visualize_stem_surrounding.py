#!/usr/bin/python

import sys, os
import random as rand
from optparse import OptionParser

import numpy as np
import scipy.stats as ss
import scipy.ndimage as sn

import corgy.builder.config as cbc
import corgy.builder.energy as cbe
import corgy.exp.kde as cek
import corgy.graph.bulge_graph as cgb
import corgy.graph.graph_pdb as cgg
import corgy.utilities.colormap as cuc
import corgy.utilities.debug as cud
import corgy.utilities.vector as cuv
import corgy.visual.pymol as cvp

import matplotlib.pyplot as plt
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
    parser.add_option('-i', '--image-filtering', dest='image_filtering', default=False, action='store_true', help='Use image filtering as a proxy for a gaussian kde.')

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
    img_sets = []
    point_sets = []
    val_sets = []
    kernels = []

    min_dims_set = []
    max_points = 80000

    for i in range(len(args)):
        stats = pa.read_csv(args[i],header=None, sep=' ')
        points = stats[['X.3', 'X.4', 'X.5']].as_matrix()

        if len(points) > max_points:
            points = np.array(rand.sample(points, max_points))

        point_sets += [points]

        if options.image_filtering:
            res = 2.

            cud.pv('points.shape')
            min_dims = np.array([min(points[:,j]) for j in xrange(points.shape[1])])
            max_dims = np.array([max(points[:,j]) for j in xrange(points.shape[1])])

            n_points = [int((max_dims[j] - min_dims[j]) / float(res)) + 1 for j in range(points.shape[1])]

            img = np.zeros(n_points)
            for p in points:
                ixs = [int((p[j] - min_dims[j]) / res) for j in xrange(points.shape[1])]
                img[ixs[0],ixs[1],ixs[2]] += 1
            img = sn.gaussian_filter(img, (2,2,2))
            img_sets += [img]

            min_dims_set += [min_dims]

        else:
            kernel = cek.gaussian_kde(points.T)
            kernels += [kernel]

            vals = cbe.my_log(kernel(points.T))
            val_sets += [vals]

    if len(point_sets) == 2:
        if options.image_filtering:
            vals = []
            points = point_sets[0]
            for p in points:
                ixs[0] = [int((p[j] - min_dims_set[0][j]) / res) for j in xrange(points.shape[1])]
                ixs[1] = [int((p[j] - min_dims_set[1][j]) / res) for j in xrange(points.shape[1])]
                #cud.pv('ixs')
                try:
                    vals += [cbe.my_log(img_sets[0][ixs[0][0], ixs[0][1], ixs[0][2]]) - cbe.my_log(img_sets[1][ixs[1][0], ixs[1][1], ixs[1][2] ])]
                except IndexError:
                    vals += [350.]
        else:
            vals = cbe.my_log(kernel_sets[0](point_sets[0])) - cbe.my_log(kernel_sets[1](point_sets[0]))
    else:
        if options.image_filtering:
            vals = []
            points = point_sets[0]
            for p in points:
                ixs = [int((p[j] - min_dims[j]) / res) for j in xrange(points.shape[1])]
                vals += [np.log(img_sets[0][ixs[0], ixs[1], ixs[2]])]
        else:
            vals = np.log(kernel_sets[0](point_sets[0]))

    mi = ss.scoreatpercentile(vals, 10)
    mx = ss.scoreatpercentile(vals, 90)
    cud.pv('mi')
    cud.pv('mx')
    #cud.pv('vals')

    #plt.hist(vals)
    #plt.show()

    for k in range(1):
        for j in range(len(point_sets[0])):
            #print "p:", p, "kernel(p):", kernel(p)
            val = vals[j]
            pp.add_sphere(point_sets[k][j], colors[i], 0.3, '', color_rgb=cuc.floatRgb(val, mi, mx))

    pp.output_pymol_file()

if __name__ == '__main__':
    main()

