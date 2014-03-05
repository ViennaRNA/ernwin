#!/usr/bin/python

from optparse import OptionParser
from bisect import bisect
import copy
import os.path as op

import forgi.utilities.debug as cud
import fess.builder.sampling as cbs
import forgi.threedee.model.stats as ftms
import fess.builder.config as cbc

import forgi.threedee.utilities.graph_pdb as cgg

from random import sample, random, seed, randint
from numpy import allclose, seterr

import forgi.threedee.model.coarse_grain as ftmc
import fess.builder.models as fbm
import forgi.threedee.utilities.rmsd as ftur

from fess.builder.sampling import StatisticsPlotter, GibbsBGSampler, SamplingStatistics
from forgi.threedee.utilities.vector import get_vector_centroid, center_on_centroid

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
        cud.pv('lines')
        for line in lines:
            if line.strip() == '':
                continue

            if line[0] == '>':
                bg = ftmc.CoarseGrainRNA()
                bg.name = line[1:].strip()
                counter = 0
            if counter % 3 == 1:
                bg.seq = line.strip()
                bg.seq_length = len(bg.seq)
            if counter % 3 == 2:
                seq = bg.seq
                name = bg.name
                bg.from_dotbracket(line.strip(), dissolve_length_one_stems=True)
                bg.seq = seq
                bg.name = name
                bgs += [bg]

            counter += 1

            cud.pv('bg.seq')
    return  bgs

def predict(bg, energies_to_sample, options):
    sm = fbm.SpatialModel(bg)

    if options.cheating:
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.CheatingEnergy(sm.bg)])]

    if not os.path.exists(options.output_dir):
        os.makedirs(options.output_dir)
    if options.output_file == None or options.output_file == sys.stdout:
        options.output_file = sys.stdout
    else:
        options.output_file = open(options.output_file, 'w')


    cbc.Configuration.sampling_output_dir = op.join(options.output_dir, bg.name)

    if options.output_dir_suffix != None:
        cbc.Configuration.sampling_output_dir = op.join(cbc.Configuration.sampling_output_dir, options.output_dir_suffix)

    if not os.path.exists(cbc.Configuration.sampling_output_dir):
        os.makedirs(cbc.Configuration.sampling_output_dir)

    if options.log_to_file:
        options.output_file = open(op.join(cbc.Configuration.sampling_output_dir, 'log.txt'), 'w')

    if options.eval_energy:
        for s in sm.bg.stem_iterator():
            cgg.add_virtual_residues(sm.bg, s)

        for energy in energies_to_sample:
            cud.pv('energy.eval_energy(sm, verbose=True)')
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
        stat = SamplingStatistics(sm, plotter, color, silent=silent, output_file=options.output_file, save_n_best = options.save_n_best, dist1 = options.dist1, dist2 = options.dist2)
        stat.step_save = options.step_save

        if options.mcmc_sampler:
            sm = fbm.SpatialModel(copy.deepcopy(bg))
            sm.constraint_energy = cbe.CombinedEnergy([cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy()])
            sm.junction_constraint_energy = cbe.RoughJunctionClosureEnergy()
            #sm.constraint_energy = cbe.CombinedEnergy([cbe.RoughJunctionClosureEnergy()])
            #sm.constraint_energy = cbe.CombinedEnergy([cbe.StemVirtualResClashEnergy()])
            #sm.constraint_energy = cbe.CombinedEnergy([cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy()])

            samplers += [cbs.MCMCSampler(sm, energy, stat, options.stats_type, options.no_rmsd)]
        else:
            sm = fbm.SpatialModel(copy.deepcopy(bg))
            sm.constraint_energy = cbe.StemVirtualResClashEnergy()
            samplers += [GibbsBGSampler(sm, energy, stat)]
        silent = True

    cud.pv('samplers')
    for i in range(options.iterations):
        if options.single_sampler:
            samplers[0].step()
        else:
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
    parser.add_option('', '--rog', dest='radius_of_gyration', default=False, action='store_true', help='Use the radius of gyration energy')
    parser.add_option('', '--encompassing-cylinder-loop-rog', dest='encompassing_cylinder_loop_radius_of_gyration', default=False, action='store_true', help='Use the cylinder_intersection and radius of gyration energy')
    parser.add_option('', '--encompassing-cylinder-rog', dest='encompassing_cylinder_radius_of_gyration', default=False, action='store_true', help='Use the cylinder_intersection and radius of gyration energy')
    parser.add_option('', '--cylinder-rog', dest='cylinder_radius_of_gyration', default=False, action='store_true', help='Use the cylinder_intersection and radius of gyration energy')
    parser.add_option('', '--cylinder-loop-rog', dest='cylinder_loop_radius_of_gyration', default=False, action='store_true', help='Use the radius of gyration energy')
    parser.add_option('', '--mloop-iloop-cylinder-loop-rog', dest='mloop_iloop_cylinder_loop_radius_of_gyration', default=False, action='store_true', help='Use the multiloop radius of gyration energy.')
    parser.add_option('', '--iloop-cylinder-loop-rog', dest='iloop_cylinder_loop_radius_of_gyration', default=False, action='store_true', help='Use the interior loop radius of gyration energy.')
    parser.add_option('', '--iloop-loop-rog', dest='iloop_loop_radius_of_gyration', default=False, action='store_true', help='Use the interior loop radius of gyration energy.')
    parser.add_option('', '--coaxiality-cylinder-loop-rog', dest='coaxiality_cylinder_loop_radius_of_gyration', default=False, action='store_true', help='Use the coaxiality, loop, and radius of gyration energy')
    parser.add_option('', '--pairwise-coaxiality-cylinder-loop-rog', dest='pairwise_coaxiality_cylinder_loop_radius_of_gyration', default=False, action='store_true', help='Use the pairwise coaxiality, loop, and radius of gyration energy')
    parser.add_option('', '--loop-rog', dest='loop_radius_of_gyration', default=False, action='store_true', help='Use the radius of gyration energy')
    parser.add_option('', '--constant-energy', dest='constant_energy', default=False, action='store_true', help='Use a constant energy')
    parser.add_option('', '--random-energy', dest='random_energy', default=False, action='store_true', help='Use a random energy')
    parser.add_option('', '--cylinder-loop', dest='cylinder_loop', default=False, action='store_true', help='Use the radius of gyration energy')
    parser.add_option('', '--stem-loop-rog', dest='stem_loop_radius_of_gyration', default=False, action='store_true', help='Use the radius of gyration energy')
    parser.add_option('', '--simple-rog', dest='simple_radius_of_gyration', default=False, action='store_true', help='Use the simple radius of gyration energy')
    parser.add_option('', '--radius-of-gyration1', dest='radius_of_gyration1', default=False, action='store_true', help='Use the radius of gyration energy')
    parser.add_option('', '--stem-stem012', dest='stem_stem012', default=False, action='store_true', help='Use the stem-stem orientation energy')
    parser.add_option('-y', '--cylinder-intersection', dest='cyl_intersect', default=False, action='store_true', help='Use the cylinder-intersection energy')
    parser.add_option('-g', '--cheating', dest='cheating', default=False, action='store_true', help='Use the rmsd from the real structure as the energy.')
    parser.add_option('', '--sequence-file', dest='sequence_file', default='', help='The file containing sequence for the structure. To be used with the --secondary-structure flag', type='str')
    parser.add_option('', '--sequence-str', dest='sequence_str', default='', help='The sequence of the structure. To be used with the --secondary-structure flag', type='str')
    parser.add_option('', '--eval-energy', dest='eval_energy', default=False, action='store_true', help='Evaluate the energy of the parameter')
    parser.add_option('', '--output-dir', dest='output_dir', default='.', help='Directory to store the sampled_structures', type='str')
    parser.add_option('', '--output-file', dest='output_file', default=None, help='File to output the information about the sampling to. Defaults to standard out', type=str)
    parser.add_option('', '--log-to-file', dest='log_to_file', default=False, help='Print a log of the output to a file in the directory where the best structures are stored.', action="store_true")

    parser.add_option('', '--save-n-best', dest='save_n_best', default=3, help='Save the best n structures.', type=int)
    parser.add_option('', '--step-save', dest='step_save', default=0, help="Save the structure at every n'th step.", type='int')
    parser.add_option('', '--loop-energy', dest='loop_energy', default=False, action='store_true', help="Add an energy function for the loop-loop interactions")
    parser.add_option('', '--loop-stem-energy', dest='loop_stem_energy', default=False, action='store_true', help="Add an energy function for the loop-loop interactions")
    parser.add_option('', '--n-loop-energy', dest='n_loop_energy', default=False, action='store_true', help="Add an energy function for the loop-loop interactions")
    parser.add_option('', '--ipe', dest='ipe', default=False, action='store_true', help="Use the interaction probability energy.")
    parser.add_option('', '--sipe', dest='sipe', default=False, action='store_true', help="Use the interaction probability energy.")
    parser.add_option('', '--no-background', dest='background', default=True, action='store_false', help="Don't use the background probability distribution.")
    parser.add_option('', '--stem-stem0-data', dest='stem_stem0_data', help='The location of the sampled stem-stem0 data', type='str', default='~/projects/ernwin/fess/stats/stem_stem_orientations_sampled.csv')
    parser.add_option('', '--stats-file', dest='stats_file', 
                      default=cbc.Configuration.stats_file, help='Use a different set of statistics for sampling', type='str') 
    parser.add_option('', '--output-dir-suffix', dest='output_dir_suffix', default=None, help="Specify an addition to the output directory", type='str')
    parser.add_option('', '--stats-type', dest='stats_type', default=None, help="Use these types of statistics.", type='str')

    parser.add_option('', '--single-sampler', dest='single_sampler', 
                      default=False, help='Use only a single sampler', action='store_true')
    parser.add_option('', '--no-rmsd', dest='no_rmsd', 
                      default=False, help='Refrain from trying to calculate the rmsd.', action='store_true')
    parser.add_option('', '--dist1', dest='dist1', default=None, help="Calculate the distance between this residue and the residue at position dist2 at every iteration", type='int')
    parser.add_option('', '--dist2', dest='dist2', default=None, help="Calculate the distance between this residue and the residue at position dist1 at every iteration", type='int')
    (options, args) = parser.parse_args()

    cud.pv('options.stem_loop_radius_of_gyration')


    if len(args) < 1:
        print "Usage: ./gibbs.py temp.comp"
        print "Or ./gibb.py temp.fa. If the extension of the argument file ends in .fa, then treat it as a fasta file."

        sys.exit(1)

    cud.pv('args')

    if options.stats_file != '':
        stats = ftms.get_angle_stats(options.stats_file)
        stats = ftms.get_loop_stats(options.stats_file)
        stats = ftms.get_stem_stats(options.stats_file)
        cbc.Configuration.stats_file = options.stats_file
        print >>sys.stderr, "options.stats_file:", cbc.Configuration.stats_file

    energies_to_sample = []
    if options.n_loop_energy:
        #energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.NLoopLoopEnergy(), cbe.StemStemOrientationEnergy([2])])]
        #energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.NLoopLoopEnergy(), cbe.NLoopStemEnergy()])]
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.NLoopLoopEnergy(), cbe.StemStemOrientationEnergy([2])])]
    if options.loop_energy:
        print "Using loop_energy"
        #energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.LoopLoopEnergy(), cbe.LoopJunctionEnergy(), cbe.LoopBulgeEnergy(), cbe.StemStemOrientationEnergy([2])])]
        #energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.LoopLoopEnergy(), cbe.StemStemOrientationEnergy([2])])]
        #energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.LoopLoopEnergy()])]
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.LoopLoopEnergy()])]
    if options.ipe:
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.InteractionProbEnergy()])]
    if options.sipe:
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.InteractionProbEnergy(), cbe.StemStemOrientationEnergy([2])])]
    if options.loop_stem_energy:
        print "Using loop_energy"
        #energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.LoopLoopEnergy(), cbe.LoopJunctionEnergy(), cbe.LoopBulgeEnergy(), cbe.StemStemOrientationEnergy([2])])]
        #energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.LoopLoopEnergy(), cbe.StemStemOrientationEnergy([2])])]
        #energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.LoopLoopEnergy()])]
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.StemStemOrientationEnergy([2])])]
    if options.cyl_intersect:
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.CylinderIntersectionEnergy()])]
    if options.stem_stem:
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.StemStemOrientationEnergy([0])])]
    if options.stem_stem0:
        sse = cbe.StemStemOrientationEnergy([0])
        sse.fake_data_location = op.expanduser(options.stem_stem0_data)
        sse.max_dist = 1000.
        sse.max_lateral_dist = 1000.
        print >>sys.stderr, 'sse'
        energies_to_sample += [cbe.CombinedEnergy([], [sse])]
    if options.radius_of_gyration:
        sse = cbe.RadiusOfGyrationEnergy()
        sse.background = options.background
        energies_to_sample += [cbe.CombinedEnergy([], [sse])]

    if options.constant_energy:
        ce = cbe.ConstantEnergy()
        energies_to_sample += [cbe.CombinedEnergy([], [ce])]

    if options.random_energy:
        re = cbe.RandomEnergy()
        energies_to_sample += [cbe.CombinedEnergy([], [re])]

    if options.pairwise_coaxiality_cylinder_loop_radius_of_gyration:
        print >>sys.stderr, "Using the coaxiality energy."
        cae = cbe.PairwiseCoaxialityEnergy()
        cie = cbe.CylinderIntersectionEnergy()
        lle = cbe.LoopLoopEnergy()
        rog = cbe.RadiusOfGyrationEnergy()
        rog.background = options.background
        energies_to_sample += [cbe.CombinedEnergy([], [lle, rog, cie, cae])]


    if options.coaxiality_cylinder_loop_radius_of_gyration:
        print >>sys.stderr, "Using the coaxiality energy."
        cae = cbe.CoaxialityEnergy()
        cie = cbe.CylinderIntersectionEnergy()
        lle = cbe.LoopLoopEnergy()
        rog = cbe.RadiusOfGyrationEnergy()
        rog.background = options.background
        energies_to_sample += [cbe.CombinedEnergy([], [lle, rog, cie, cae])]

    if options.iloop_loop_radius_of_gyration:
        lle = cbe.LoopLoopEnergy()
        ile = cbe.LoopILoopEnergy()
        rog = cbe.RadiusOfGyrationEnergy()
        rog.background = options.background
        energies_to_sample += [cbe.CombinedEnergy([], [lle, rog, ile])]

    if options.mloop_iloop_cylinder_loop_radius_of_gyration:
        mle = cbe.LoopMLoopEnergy()
        cie = cbe.CylinderIntersectionEnergy()
        lle = cbe.LoopLoopEnergy()
        ile = cbe.LoopILoopEnergy()
        rog = cbe.RadiusOfGyrationEnergy()
        rog.background = options.background
        energies_to_sample += [cbe.CombinedEnergy([], [lle, rog, cie, ile, mle])]

    if options.iloop_cylinder_loop_radius_of_gyration:
        cie = cbe.CylinderIntersectionEnergy()
        lle = cbe.LoopLoopEnergy()
        ile = cbe.LoopILoopEnergy()
        rog = cbe.RadiusOfGyrationEnergy()
        rog.background = options.background
        energies_to_sample += [cbe.CombinedEnergy([], [lle, rog, cie, ile])]

    if options.cylinder_loop_radius_of_gyration:
        cie = cbe.CylinderIntersectionEnergy()
        lle = cbe.LoopLoopEnergy()
        rog = cbe.RadiusOfGyrationEnergy()
        rog.background = options.background
        energies_to_sample += [cbe.CombinedEnergy([], [lle, rog, cie])]

    if options.cylinder_loop:
        lle = cbe.LoopLoopEnergy()
        cie = cbe.CylinderIntersectionEnergy()
        sse.background = options.background
        energies_to_sample += [cbe.CombinedEnergy([], [cie, lle])]

    if options.encompassing_cylinder_loop_radius_of_gyration:
        ece = cbe.EncompassingCylinderEnergy()
        sse = cbe.RadiusOfGyrationEnergy()
        lle = cbe.LoopLoopEnergy()
        cie = cbe.CylinderIntersectionEnergy()
        sse.background = options.background
        energies_to_sample += [cbe.CombinedEnergy([], [sse, cie, ece, lle])]

    if options.encompassing_cylinder_radius_of_gyration:
        ece = cbe.EncompassingCylinderEnergy()
        sse = cbe.RadiusOfGyrationEnergy()
        cie = cbe.CylinderIntersectionEnergy()
        sse.background = options.background
        energies_to_sample += [cbe.CombinedEnergy([], [sse, cie, ece])]

    if options.cylinder_radius_of_gyration:
        sse = cbe.RadiusOfGyrationEnergy()
        cie = cbe.CylinderIntersectionEnergy()
        sse.background = options.background
        energies_to_sample += [cbe.CombinedEnergy([], [sse, cbe.CylinderIntersectionEnergy()])]

    if options.loop_radius_of_gyration:
        sse = cbe.RadiusOfGyrationEnergy()
        sse.background = options.background
        energies_to_sample += [cbe.CombinedEnergy([], [sse, cbe.LoopLoopEnergy()])]

    if options.stem_loop_radius_of_gyration:
        sse1 = cbe.StemStemOrientationEnergy([2])
        sse1.max_dist = 40
        sse1.max_lateral_dist = 12.

        sse1.real_data = sse1.load_stem_stem_data('fess/stats/stem_stem_orientations.csv')
        sse1.fake_data = sse1.load_stem_stem_data('fess/stats/stem_stem_orientations_random_loop_radius_gyration_beta_29.csv')
        sse = cbe.RadiusOfGyrationEnergy()
        sse.background = options.background
        energies_to_sample += [cbe.CombinedEnergy([], [sse, cbe.LoopLoopEnergy(), sse1])]


    if options.simple_radius_of_gyration:
        sse = cbe.RadiusOfGyrationEnergy(dist_type="beta", adjustment=0.6)
        sse.background = options.background
        energies_to_sample += [cbe.CombinedEnergy([], [sse])]
        
    if options.radius_of_gyration1:
        sse = cbe.RadiusOfGyrationEnergy()
        sse.background = False
        sse1 = cbe.RadiusOfGyrationEnergy()
        sse1.background = True
        sse1.sampled_stats_fn = 'fess/stats/subgraph_radius_of_gyration_target.csv'
        energies_to_sample += [cbe.CombinedEnergy([], [sse, sse1])]

    if options.stem_stem2:
        #energies_to_sample += [cbe.CombinedEnergy([], [cbe.CoarseStemClashEnergy(), cbe.StemVirtualResClashEnergy(), cbe.RoughJunctionClosureEnergy(), cbe.StemStemOrientationEnergy([2]), cbe.CylinderIntersectionEnergy()])]
        sse0 = cbe.StemStemOrientationEnergy([0])
        sse0.max_dist = 1000.
        sse0.max_lateral_dist = 1000.
        sse0.beta = True

        sse2 = cbe.StemStemOrientationEnergy([2])

        print "Using stem_stem2 energy"
        energies_to_sample += [cbe.CombinedEnergy([], [sse2])]

    if options.stem_stem02:
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.StemStemOrientationEnergy([0]),cbe.StemStemOrientationEnergy([2])])]
    if options.stem_stem012:
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.StemStemOrientationEnergy([0,1,2])])]
    if options.helix_orientation:
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.ImgHelixOrientationEnergy()])]
        #energies_to_sample += [cbe.CombinedEnergy([], [cbe.StemVirtualResClashEnergy(), cbe.ImgHelixOrientationEnergy()])]
    if options.constrained_energy:
        energies_to_sample += [cbe.CombinedEnergy([], [cbe.RandomEnergy()])]
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

    bgs = []

    for arg in args:
        if arg[-3:] == '.fa':
            bgs += bgs_from_fasta(arg)
        else:
            bgs += [ftmc.CoarseGrainRNA(arg)]

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

