#!/usr/bin/python

from optparse import OptionParser
from bisect import bisect
from copy import deepcopy

from bobbins_config import ConstructionConfig
from random import sample, random, seed
from numpy import allclose, seterr

from corgy.graph.bulge_graph import BulgeGraph
from corgy.builder.energy import LongRangeInteractionCount, RandomEnergy
from corgy.builder.energy import CombinedEnergy
from corgy.builder.models import SpatialModel
from corgy.builder.rmsd import centered_rmsd

from corgy.builder.sampling import StatisticsPlotter, GibbsBGSampler, SamplingStatistics

from corgy.utilities.vector import get_vector_centroid, center_on_centroid

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
    parser.add_option('-a', '--angle_stats', dest='angle_stats_fn', default=ConstructionConfig.angle_stats_file, help='Location of the angle statistics file.') 
    parser.add_option('-p', '--plot', dest='plot', default=False, action='store_true', help="Plot the energies as they are encountered")

    (options, args) = parser.parse_args()

    if len(args) < 1:
        print "Usage: ./gibbs.py temp.comp"
        sys.exit(1)

    bg = BulgeGraph(args[0])
    bg.calc_bp_distances()
    sm = SpatialModel(bg)

    '''
    ef1 = pickle.load(open('energies/%s/LongRangeInteractionCount' % (bg.name), 'r'))
    ef2 = pickle.load(open('energies/%s/LongRangeDistanceEnergy' % (bg.name), 'r'))
    ef3 = pickle.load(open('energies/%s/JunctionClosureEnergy' % (bg.name), 'r'))
    ef4 = pickle.load(open('energies/%s/SkewNormalInteractionEnergy' % (bg.name), 'r'))
    '''

    #print "LRIC.target", 

    #energy_function = CombinedEnergy([ef1, ef4, ef3])
    #energy_function = CombinedEnergy([ef1, ef3, ef4])
    #energy_function = CombinedEnergy([ef4])
    #energy_function = ef3

    energy_function = pickle.load(open('/home/mescalin/pkerp/projects/ernwin/bobbins/energy/%s/1000/SkewNormalInteractionEnergy/LongRangeInteractionCount/JunctionClosureEnergy/CombinedEnergy.energy' % (bg.name), 'r'))

    if options.plot:
        plotter = StatisticsPlotter()
    else:
        plotter = None

    #plotter = None

    stats = SamplingStatistics(sm, plotter, 'b', silent=False)
    random_stats = SamplingStatistics(sm, plotter, 'r', silent=True)

    gs = GibbsBGSampler(deepcopy(sm), energy_function, stats)
    gs_random = GibbsBGSampler(deepcopy(sm), CombinedEnergy([RandomEnergy()]), random_stats)

    for i in range(options.iterations):
        gs.step()
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
        pdb.post_mortem()
        raise

