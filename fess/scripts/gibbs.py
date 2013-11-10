#!/usr/bin/python

from optparse import OptionParser
from bisect import bisect
import copy
import os.path as op

import forgi.utilities.debug as cud
import fess.builder.sampling as cbs
import fess.builder.config as cbc

import forgi.threedee.utilities.graph_pdb as cgg

from random import sample, random, seed, randint
from numpy import allclose, seterr

from forgi.graph.bulge_graph import BulgeGraph
from fess.builder.models import SpatialModel
from forgi.threedee.utilities.rmsd import centered_rmsd

from fess.builder.sampling import StatisticsPlotter, GibbsBGSampler, SamplingStatistics
from forgi.threed.utilities.vector import get_vector_centroid, center_on_centroid

import fess.builder.config as conf
import fess.builder.energy as cbe
import os

import sys
import pickle, pdb

from math import exp

from sys import stderr

def draw_helper():
    draw()
    pass

def bgs_from_fasta(fasta_file):
    bgs = []

    if not op.exists(fasta_file):
        print >>sys.stderr, "The specified fasta file does not exist: %s" \
                % (fasta_file)

    with open(fasta_file, 'r') as f:
        # assume there is only one sequence and dotplot in the fastdp file
        lines = f.readlines()
        for line in lines:
            if line.strip() == '':
                continue

            if line[0] == '>':
                bg = BulgeGraph()
                bg.name = line[1:].strip()
                counter = 0
            if counter % 3 == 1:
                bg.seq = line.strip()
                bg.length = len(bg.seq)
            if counter % 3 == 2:
                bg.from_dotbracket(line.strip())
                bgs += [bg]

            counter += 1
    return  bgs

def predict(bg, energies_to_sample, options):
    sm = SpatialModel(bg)

    if not os.path.exists(options.output_dir):
        os.makedirs(options.output_dir)
    if options.output_file == None or options.output_file == sys.stdout:
        options.output_file = sys.stdout
    else:
        options.output_file = open(options.output_file, 'w')

    cbc.Configuration.sampling_output_dir = op.join(options.output_dir, bg.name)
    if not os.path.exists(cbc.Configuration.sampling_output_dir):
        os.makedirs(cbc.Configuration.sampling_output_dir)

    if options.eval_energy:
        for s in sm.bg.stems():
            cgg.add_virtual_residues(sm.bg, s)

        for energy in energies_to_sample:
            cud.pv('energy.eval_energy(sm)')
        sys.exit(1)

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
        stat = SamplingStatistics(sm, plotter, color, silent=silent, output_file=options.output_file, save_n_best = options.save_n_best)
        stat.step_save = options.step_save

        if options.mcmc_sampler:
            sm = SpatialModel(copy.deepcopy(bg))
            sm.constraint_energy = cbe.CombinedEnergy([cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy()])
            sm.junction_constraint_energy = cbe.RoughJunctionClosureEnergy()
            #sm.constraint_energy = cbe.CombinedEnergy([cbe.RoughJunctionClosureEnergy()])
            #sm.constraint_energy = cbe.CombinedEnergy([cbe.StemVirtualResClashEnergy()])
            #sm.constraint_energy = cbe.CombinedEnergy([cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy()])
            samplers += [cbs.MCMCSampler(sm, energy, stat)]
        else:
            sm = SpatialModel(copy.deepcopy(bg))
            sm.constraint_energy = cbe.StemVirtualResClashEnergy()
            samplers += [GibbsBGSampler(sm, energy, stat)]
        silent = True

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
def main():
    #seed(2)
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
    parser.add_option('', '--stem-stem0', dest='stem_stem0', default=False, action='store_true', help='Use the stem-stem orientation energy')
    parser.add_option('', '--stem-stem2', dest='stem_stem2', default=False, action='store_true', help='Use the stem-stem orientation energy')
    parser.add_option('', '--stem-stem02', dest='stem_stem02', default=False, action='store_true', help='Use the stem-stem orientation energy')
    parser.add_option('', '--stem-stem012', dest='stem_stem012', default=False, action='store_true', help='Use the stem-stem orientation energy')
    parser.add_option('-y', '--cylinder-intersection', dest='cyl_intersect', default=False, action='store_true', help='Use the cylinder-intersection energy')
    parser.add_option('-g', '--cheating', dest='cheating', default=False, action='store_true', help='Use the rmsd from the real structure as the energy.')
    parser.add_option('', '--secondary-structure', dest='secondary_structure', default=False, action='store_true', help='Take a secondary structure as input instead of a bulge graph')
    parser.add_option('', '--sequence-file', dest='sequence_file', default='', help='The file containing sequence for the structure. To be used with the --secondary-structure flag', type='str')
    parser.add_option('', '--sequence-str', dest='sequence_str', default='', help='The sequence of the structure. To be used with the --secondary-structure flag', type='str')
    parser.add_option('', '--eval-energy', dest='eval_energy', default=False, action='store_true', help='Evaluate the energy of the parameter')
    parser.add_option('', '--output-dir', dest='output_dir', default='.', help='Directory to store the sampled_structures', type='str')
    parser.add_option('', '--output-file', dest='output_file', default=None, help='File to output the information about the sampling to. Defaults to standard out', type=str)
    parser.add_option('', '--save-n-best', dest='save_n_best', default=3, help='Save the best n structures.', type=int)
    parser.add_option('', '--step-save', dest='step_save', default=False, action='store_true', help="Save the structure at each step.")
    parser.add_option('', '--loop-energy', dest='loop_energy', default=False, action='store_true', help="Add an energy function for the loop-loop interactions")
    parser.add_option('', '--loop-stem-energy', dest='loop_stem_energy', default=False, action='store_true', help="Add an energy function for the loop-loop interactions")
    parser.add_option('', '--n-loop-energy', dest='n_loop_energy', default=False, action='store_true', help="Add an energy function for the loop-loop interactions")
    parser.add_option('', '--ipe', dest='ipe', default=False, action='store_true', help="Use the interaction probability energy.")
    parser.add_option('', '--sipe', dest='sipe', default=False, action='store_true', help="Use the interaction probability energy.")
    parser.add_option('', '--fasta', dest='fasta', default='', help="Specify a fastdb file containing an identifier, a sequence and a dotbracket string indicating the secondary structure.", type='str')
    parser.add_option('', '--stem-stem0-data', dest='stem_stem0_data', help='The location of the sampled stem-stem0 data', type='str', default='~/projects/ernwin/fess/stats/stem_stem_orientations_sampled.csv')

    (options, args) = parser.parse_args()

    if len(args) < 1 and options.secondary_structure is False and options.fasta == '':
        print "Usage: ./gibbs.py temp.comp"
        sys.exit(1)

    cud.pv('args')
    cud.pv('options.secondary_structure')

    energies_to_sample = []
    if options.n_loop_energy:
        #energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.NLoopLoopEnergy(), cbe.StemStemOrientationEnergy([2])])]
        #energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.NLoopLoopEnergy(), cbe.NLoopStemEnergy()])]
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.NLoopLoopEnergy(), cbe.StemStemOrientationEnergy([2])])]
    if options.loop_energy:
        print "Using loop_energy"
        #energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.LoopLoopEnergy(), cbe.LoopJunctionEnergy(), cbe.LoopBulgeEnergy(), cbe.StemStemOrientationEnergy([2])])]
        #energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.LoopLoopEnergy(), cbe.StemStemOrientationEnergy([2])])]
        #energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.LoopLoopEnergy()])]
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.LoopLoopEnergy()])]
    if options.ipe:
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.InteractionProbEnergy()])]
    if options.sipe:
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.InteractionProbEnergy(), cbe.StemStemOrientationEnergy([2])])]
    if options.loop_stem_energy:
        print "Using loop_energy"
        #energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.LoopLoopEnergy(), cbe.LoopJunctionEnergy(), cbe.LoopBulgeEnergy(), cbe.StemStemOrientationEnergy([2])])]
        #energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.LoopLoopEnergy(), cbe.StemStemOrientationEnergy([2])])]
        #energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.LoopLoopEnergy()])]
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.LoopLoopEnergy(), cbe.StemStemOrientationEnergy([2])])]
    if options.cheating:
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.CheatingEnergy(sm.bg)])]
    if options.cyl_intersect:
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.CylinderIntersectionEnergy()])]
    if options.stem_stem:
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.StemStemOrientationEnergy([0])])]
    if options.stem_stem0:
        sse = cbe.StemStemOrientationEnergy([0])
        sse.fake_data_location = op.expanduser(options.stem_stem0_data)
        sse.max_dist = 1000.
        sse.max_lateral_dist = 1000.
        print >>sys.stderr, 'sse'
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), sse])]
    if options.stem_stem2:
        #energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.StemStemOrientationEnergy([2]), cbe.CylinderIntersectionEnergy()])]
        sse0 = cbe.StemStemOrientationEnergy([0])
        sse0.max_dist = 1000.
        sse0.max_lateral_dist = 1000.
        sse0.beta = True


        sse2 = cbe.StemStemOrientationEnergy([2])

        print "Using stem_stem2 energy"
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), sse2])]

    if options.stem_stem02:
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.StemStemOrientationEnergy([0]),cbe.StemStemOrientationEnergy([2])])]
    if options.stem_stem012:
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.StemStemOrientationEnergy([0,1,2])])]
    if options.helix_orientation:
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.ImgHelixOrientationEnergy()])]
        #energies_to_sample += [cbe.CombinedEnergy([], [cbe.StemVirtualResClashEnergy(), cbe.ImgHelixOrientationEnergy()])]
    if options.constrained_energy:
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.RandomEnergy(), cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy()])]
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

    if options.fasta != '':
        bgs = bgs_from_fasta(options.fasta)

    elif options.secondary_structure:
        if options.sequence_file == '' and options.sequence_str != '':
            print >>sys.stderr, "Sequence needs to be provided with --sequence"
        print >>sys.stderr, "Secondary structure provided in lieu of a bulge-graph"
        bg = BulgeGraph()

        if options.sequence_file != '':
            with open(options.sequence, 'r') as f:
                seq = "".join(f.readlines())
                bg.seq = seq.upper()
        else:
            bg.seq = options.sequence_str

        bg.from_dotbracket_file(args[0])
        bgs = [bg]
    else:
        bgs = [BulgeGraph(args[0])]

    #bg.calc_bp_distances()
    for bg in bgs:
        predict(bg, energies_to_sample, options)

if __name__ == '__main__':
    seed_num = randint(0, sys.maxint)
    try:
        seed(seed_num)
        main()
    except:
        #pdb.post_mortem()
        print >>sys.stderr, "random seed:", seed_num
        raise

