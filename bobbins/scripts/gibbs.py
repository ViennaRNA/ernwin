#!/usr/bin/python

from optparse import OptionParser
from bisect import bisect
import copy

from random import sample, random, seed
from numpy import allclose, seterr



from corgy.graph.bulge_graph import BulgeGraph
from corgy.builder.models import SpatialModel
from corgy.builder.rmsd import centered_rmsd

from corgy.builder.sampling import StatisticsPlotter, GibbsBGSampler, SamplingStatistics

from corgy.utilities.vector import get_vector_centroid, center_on_centroid

import corgy.builder.config as conf
import corgy.builder.energy as cbe
import os

import sys
import pickle, pdb

from math import exp

from sys import stderr
from pylab import plot,show, hist, xlabel, ylabel, xlim, ylim, ion, draw, ioff, clf

def draw_helper():
    draw()
    pass

def main():
    seed(2)
    ion()
    #seterr(all='ignore')
    #seterr(all='raise')
    parser = OptionParser()

    parser.add_option('-e', '--energy', dest='energy', default='energies/lrde.energy', help="The energy function to use when evaluating structures")
    parser.add_option('-i', '--iterations', dest='iterations', default=10, help='Number of structures to generate', type='int')
    parser.add_option('-b', '--best_filename', dest='best_filename', default='best.coord', help="The filename to dump the best (least rmsd structure) into", type='str')
    parser.add_option('-p', '--plot', dest='plot', default=False, action='store_true', help="Plot the energies as they are encountered")
    parser.add_option('-d', '--distance', dest='distance_energy', default=False, action='store_true', help='Use the DistanceEnergy energy')
    parser.add_option('-r', '--random', dest='random_energy', default=False, action='store_true', help='Use the RandomEnergy energy')
    parser.add_option('-s', '--step-random', dest='step_random', default=False, action='store_true', help='Concurrently sample with a random energy.')
    parser.add_option('-o', '--helix_orientation', dest='helix_orientation', default=False, action='store_true', help='Sample using the helix orientation energy')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        print "Usage: ./gibbs.py temp.comp"
        sys.exit(1)

    bg = BulgeGraph(args[0])
    sm = SpatialModel(bg)

    if options.random_energy:
        energy_function = cbe.CombinedEnergy([], [cbe.RandomEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy()])
    elif options.distance_energy:
        energy_function = cbe.DistanceEnergy(bg.get_long_range_constraints())
    elif options.helix_orientation:
        energy_function = cbe.CombinedEnergy([], [cbe.StemVirtualResClashEnergy(), cbe.ImgHelixOrientationEnergy()])
    else:
        bg.calc_bp_distances()
        energy_function = pickle.load(open(os.path.join(conf.Configuration.base_dir, 'bobbins/energy/%s/1000/SkewNormalInteractionEnergy/LongRangeInteractionCount/JunctionClosureEnergy/CombinedEnergy.energy' % (bg.name)), 'r'))



    if options.plot:
        plotter = StatisticsPlotter()
    else:
        plotter = None

    #plotter = None

    stats = SamplingStatistics(sm, plotter, 'b', silent=False)
    random_stats = SamplingStatistics(sm, plotter, 'r', silent=True)

    gs = GibbsBGSampler(sm, energy_function, stats)
    gs_random = GibbsBGSampler(copy.deepcopy(sm), cbe.CombinedEnergy([], [cbe.RandomEnergy(), cbe.StemVirtualResClashEnergy()]), random_stats)

    for i in range(options.iterations):
        gs.step()
        if options.step_random:
            gs_random.step()

    stats.print_final_stats(energy_function)
    stats.save_top()

    if plotter:
        plotter.finish()
        pass

if __name__ == '__main__':
    try:
        main()
    except:
        #pdb.post_mortem()
        raise

