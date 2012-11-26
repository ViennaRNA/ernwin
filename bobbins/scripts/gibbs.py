#!/usr/bin/python

from optparse import OptionParser
from bisect import bisect
import copy

import corgy.utilities.debug as cud
import corgy.builder.sampling as cbs

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
    parser.add_option('-c', '--constrained', dest='constrained_energy', default=False, action='store_true', help='Use the ConstrainedRandomEnergy energy')
    parser.add_option('-r', '--step-random', dest='step_random', default=False, action='store_true', help='Concurrently sample with a random energy.')
    parser.add_option('-o', '--helix_orientation', dest='helix_orientation', default=False, action='store_true', help='Sample using the helix orientation energy')
    parser.add_option('-m', '--mcmc', dest='mcmc_sampler', default=False, action='store_true', help='Sample using the mcmc sampler.')
    parser.add_option('-s', '--stem-stem', dest='stem_stem', default=False, action='store_true', help='Use the stem-stem orientation energy')
    parser.add_option('', '--secondary-structure', dest='secondary_structure', default=False, action='store_true', help='Take a secondary structure as input instead of a bulge graph')


    (options, args) = parser.parse_args()

    if len(args) < 1:
        print "Usage: ./gibbs.py temp.comp"
        sys.exit(1)

    if options.secondary_structure:
        print >>sys.stderr, "Secondary structure provided in lieu of a bulge-graph"
        bg = BulgeGraph()
        bg.from_dotbracket_file(args[0])
    else:
        bg = BulgeGraph(args[0])

    sm = SpatialModel(bg)

    energies_to_sample = []
    
    if options.stem_stem:
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.StemStemOrientationEnergy()])]
    if options.helix_orientation:
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.ImgHelixOrientationEnergy()])]
        #energies_to_sample += [cbe.CombinedEnergy([], [cbe.StemVirtualResClashEnergy(), cbe.ImgHelixOrientationEnergy()])]
    if options.constrained_energy:
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.RandomEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy()])]
        #energies_to_sample += [cbe.CombinedEnergy([], [cbe.RandomEnergy(), cbe.StemVirtualResClashEnergy()])]
    if options.distance_energy:
        energies_to_sample += [cbe.DistanceEnergy(bg.get_long_range_constraints())]
    if options.step_random:
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.RandomEnergy()])]
    '''
    else:
        bg.calc_bp_distances()
        energy_function = pickle.load(open(os.path.join(conf.Configuration.base_dir, 'bobbins/energy/%s/1000/SkewNormalInteractionEnergy/LongRangeInteractionCount/JunctionClosureEnergy/CombinedEnergy.energy' % (bg.name)), 'r'))
    '''

    if options.plot:
        plotter = StatisticsPlotter()
    else:
        plotter = None

    colors = ['g','y','r']
    samplers = []
    stats = []

    # only samples from the first energy will be saved
    silent = False

    for color,energy in zip(colors, energies_to_sample):
        stat = SamplingStatistics(sm, plotter, color, silent=silent) 
        if options.mcmc_sampler:
            samplers += [cbs.MCMCSampler(copy.deepcopy(sm), energy, stat)]
        else:
            samplers += [GibbsBGSampler(copy.deepcopy(sm), energy, stat)]
        silent = True

    cud.pv('samplers')

    for i in range(options.iterations):
        for s in samplers:
            s.step()

    #stats.print_final_stats(energy_function)
    #stats.save_top()

    if plotter:
        plotter.finish()
        #plotter.plt.ioff()
        #plt.show()
        pass

if __name__ == '__main__':
    try:
        main()
    except:
        #pdb.post_mortem()
        raise

