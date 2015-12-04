#!python

import collections as col
import copy
import random
import optparse
import os
import os.path as op
import subprocess as spr
import sys

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.stats as ftms
import forgi.threedee.utilities.graph_pdb as cgg
import forgi.utilities.debug as fud

import fess
import fess.builder.energy as fbe
import fess.builder.config as cbc
import fess.builder.models as fbm
import fess.builder.sampling as fbs

def bgs_from_fasta(fasta_file):
    bgs = []

    if not op.exists(fasta_file):
        print >>sys.stderr, "The specified fasta file does not exist: %s" \
                % (fasta_file)

    with open(fasta_file, 'r') as f:
        # assume there is only one sequence and dotplot in the fastdp file
        lines = f.read()
        bg = ftmc.CoarseGrainRNA()
        bg.from_fasta(lines, dissolve_length_one_stems=False)
        bgs += [bg]
    return  bgs

def predict(bg, energies_to_sample, options):
    fud.pv('energies_to_sample[0].energies')

    if options.cheating:
        sm = fbm.SpatialModel(bg)
        #energies_to_sample += [fbe.CombinedEnergy([], [fbe.CheatingEnergy(sm.bg)])]
        energies_to_sample = [fbe.CheatingEnergy(sm.bg)]

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

    if options.fix_all_loops:
        options.fix_loop = ','.join([d for d in bg.defines if d[0] == 'i'])

    if options.jared_dir is not None:
        # run the jar3d_annotate script to get a list of potential statistics for each interior loop
        jared_script = op.join(options.jared_dir, 'scripts/annotate_structure.py')
        jared_data = op.join(options.jared_dir, 'JAR3D')

        filtered_stats_fn = op.join(cbc.Configuration.sampling_output_dir,
                                    'filtered.stats')

        cmd = ['python', jared_script, options.bg_filename, '-o', jared_data,
               '-m', '-e', '-d', jared_data]

        fud.pv('cmd')
        p = spr.Popen(cmd, stdout=spr.PIPE)
        out, err = p.communicate()

        with open(filtered_stats_fn, 'w') as filtered_out:
            filtered_out.write(out)

        filtered_stats = ftms.FilteredConformationStats(stats_file=options.stats_file,
                                                        filter_filename=filtered_stats_fn)
        ftms.set_conformation_stats(filtered_stats)
        print >>sys.stderr, "Using JAR3D filtered stats"
    elif options.filtered_stats_file is not None:
        filtered_stats = ftms.FilteredConformationStats(stats_file=options.stats_file,
                                                        filter_filename=options.filtered_stats_file)
        ftms.set_conformation_stats(filtered_stats)
    elif options.fix_loop is not None:
        filtered_stats = ftms.FilteredConformationStats(stats_file=options.stats_file)
        filtered_stats.filtered_stats = col.defaultdict(list)

        for element_to_fix in options.fix_loop.split(','):
            print >>sys.stderr, "fixing loop", element_to_fix
            if element_to_fix[0] != 'i' and element_to_fix[0] != 'm':
                print >>sys.stderr, "Cannot fix non-interior loop or multi-loop stats, yet!"
                sys.exit(1)
            as1, as2 = bg.get_bulge_angle_stats(element_to_fix)


            filtered_stats.filtered_stats[(element_to_fix, as1.ang_type)] += [as1]
            filtered_stats.filtered_stats[(element_to_fix, as2.ang_type)] += [as2]

        ftms.set_conformation_stats(filtered_stats)
        fud.pv('element_to_fix')

    elif options.stats_file is not None:
        cbc.Configuration.stats_file = options.stats_file
        print >>sys.stderr, "1"
        ftms.set_conformation_stats(ftms.ConformationStats(options.stats_file))

    sm = fbm.SpatialModel(bg)

    if options.log_to_file:
        options.output_file = open(op.join(cbc.Configuration.sampling_output_dir, 'log.txt'), 'w')

    if options.eval_energy:
        for s in sm.bg.stem_iterator():
            cgg.add_virtual_residues(sm.bg, s)

        for energy in energies_to_sample:
            fud.pv('energy.eval_energy(sm, verbose=True, background=False)')
        sys.exit(1)

    if options.plot:
        plotter = fbs.StatisticsPlotter()
    else:
        plotter = None

    colors = ['g','y','r']
    samplers = []

    # parse the distances that we want to keep track of 
    to_track_dists = []
    if options.dists is not None:
        for distance_pair in options.dists.split(':'):
            to_track_dists += [map(int, distance_pair.split(','))]

    # only samples from the first energy will be saved
    silent = False

    for color,energy in zip(colors, energies_to_sample):
        fud.pv('options.no_rmsd')
        stat = fbs.SamplingStatistics(sm, plotter, color, silent=silent, 
                                      output_file=options.output_file, 
                                      save_n_best = options.save_n_best, 
                                      dists = to_track_dists, 
                                      save_iterative_cg_measures=options.save_iterative_cg_measures, 
                                      no_rmsd = options.no_rmsd)
        stat.step_save = options.step_save

        fud.pv('options.mcmc_sampler')
        if options.mcmc_sampler:
            sm = fbm.SpatialModel(copy.deepcopy(bg))

            sm.constraint_energy = fbe.CombinedEnergy([])
            sm.junction_constraint_energy = fbe.CombinedEnergy([])

            if not (options.cheating or options.no_constraint):
                sm.constraint_energy = fbe.CombinedEnergy([fbe.CoarseStemClashEnergy(), fbe.StemVirtualResClashEnergy()])
                sm.junction_constraint_energy = fbe.RoughJunctionClosureEnergy()
            else:
                sm.constraint_energy = None
                sm.junction_constraint_energy = None

            #sm.constraint_energy = fbe.CombinedEnergy([fbe.RoughJunctionClosureEnergy()])
            #sm.constraint_energy = fbe.CombinedEnergy([fbe.StemVirtualResClashEnergy()])
            #sm.constraint_energy = fbe.CombinedEnergy([fbe.StemVirtualResClashEnergy(), fbe.RoughJunctionClosureEnergy()])
            if options.track_energies:
                energies_to_track = [fbe.RadiusOfGyrationEnergy()]

                fud.pv('len(list(bg.hloop_iterator()))')
                if len(list(bg.hloop_iterator())) > 1:
                    energies_to_track += [fbe.ShortestLoopDistanceEnergy()]
                    for h in bg.hloop_iterator():
                        energies_to_track += [fbe.ShortestLoopDistancePerLoop(h)]

                energies_to_track += [fbe.AMinorEnergy(loop_type='h')]
                #energies_to_track += [fbe.AMinorEnergy(loop_type='i')]
            else:
                energies_to_track = []

            fud.pv('energies_to_track')
            fud.pv('energy')
            sampler = fbs.MCMCSampler(sm, energy, stat, options.stats_type, options.no_rmsd, energies_to_track=energies_to_track)
            sampler.dump_measures = options.dump_energies
            samplers += [sampler]
        else:
            sm = fbm.SpatialModel(copy.deepcopy(bg))
            sm.constraint_energy = fbe.StemVirtualResClashEnergy()
            samplers += [fbs.GibbsBGSampler(sm, energy, stat)]
        silent = True

    fud.pv('samplers')
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
    parser = optparse.OptionParser()

    parser.add_option('', '--loop-energy', dest='loop_energy', default=False, action='store_true', help='Use the radius of gyration energy')
    parser.add_option('', '--dump-energies', dest='dump_energies', default=False, action='store_true', help='Dump the energies to file')
    parser.add_option('', '--track-energies', dest='track_energies', default=False, help='Track additional energy for diagnostics', action='store_true')
    parser.add_option('', '--energy-prefactor', dest='energy_prefactor', default=30, help='A multiplier for the energy', type='int')
    parser.add_option('-e', '--energy', dest='energy', default='energies/lrde.energy', help="The energy function to use when evaluating structures")
    parser.add_option('-i', '--iterations', dest='iterations', default=10000, help='Number of structures to generate', type='int')
    parser.add_option('-b', '--best_filename', dest='best_filename', default='best.coord', help="The filename to dump the best (least rmsd structure) into", type='str')
    parser.add_option('-p', '--plot', dest='plot', default=False, action='store_true', help="Plot the energies as they are encountered")
    parser.add_option('-d', '--distance', dest='distance_energy', default=False, action='store_true', help='Use the DistanceEnergy energy')
    parser.add_option('-c', '--clamp', dest='clamp', default=None, help='Clamp two elements together (i.e. add an energy with a target distance of 10 angstroms). The energy should be formatted as p1,p2:p3,p4:p5,p6 where p1 and p2 are clamped, p3 and p4 are clamped and p5 and p6 are clamped.', type='str')
    parser.add_option('-m', '--mcmc', dest='mcmc_sampler', default=True, action='store_true', help='Sample using the mcmc sampler.')
    parser.add_option('', '--rog', dest='radius_of_gyration', default=False, action='store_true', help='Use the radius of gyration energy')
    parser.add_option('', '--cylinder-rog', dest='cylinder_radius_of_gyration', default=False, action='store_true', help='Use the cylinder_intersection and radius of gyration energy')
    parser.add_option('', '--aminor-perloop-rog', dest='aminor_perloop_radius_of_gyration', default=False, action='store_true', help='Use the aminor and radius of gyration energies')
    parser.add_option('', '--specific-aminor', dest='specific_aminor', default=None, help='Use the specific aminor energy', type='str')
    parser.add_option('', '--aminor-perloop', dest='aminor_perloop', default=False, action='store_true', help='Use the aminor and radius of gyration energies')
    parser.add_option('', '--aminor-shortestloop', dest='aminor_shortestloop', default=False, action='store_true', help='Use the aminor and radius of gyration energies')
    parser.add_option('', '--aminor-rog', dest='aminor_radius_of_gyration', default=False, action='store_true', help='Use the aminor and radius of gyration energies')
    parser.add_option('', '--aminor', dest='aminor', default=False, action='store_true', help='Use the aminor and radius of gyration energies')
    parser.add_option('', '--cylinder-perloop-rog', dest='cylinder_perloop_radius_of_gyration', default=False, action='store_true', help='Use the radius of gyration energy')
    parser.add_option('', '--cylinder-shortestloop-rog', dest='cylinder_shortestloop_radius_of_gyration', default=False, action='store_true', help='Use the radius of gyration energy')
    parser.add_option('', '--cylinder-loop-rog', dest='cylinder_loop_radius_of_gyration', default=False, action='store_true', help='Use the radius of gyration energy')
    parser.add_option('', '--loop-rog', dest='loop_radius_of_gyration', default=False, action='store_true', help='Use the radius of gyration energy')
    parser.add_option('', '--constant-energy', dest='constant_energy', default=False, action='store_true', help='Use a constant energy')
    parser.add_option('', '--random-energy', dest='random_energy', default=False, action='store_true', help='Use a random energy')
    parser.add_option('', '--cylinder-loop', dest='cylinder_loop', default=False, action='store_true', help='Use the radius of gyration energy')
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
    parser.add_option('', '--no-background', dest='background', default=True, action='store_false', help="Don't use the background probability distribution.")
    parser.add_option('', '--stats-file', dest='stats_file', 
                      default=fess.data_file('stats/combined.stats'), help='Use a different set of statistics for sampling', type='str') 
    parser.add_option('', '--filtered-stats-file', dest='filtered_stats_file', 
                      default=None, 
                      help='Filter the statistics used for sampling using some other file.', type='str') 
    parser.add_option('', '--output-dir-suffix', dest='output_dir_suffix', default=None, help="Specify an addition to the output directory", type='str')
    parser.add_option('', '--stats-type', dest='stats_type', default=None, help="Use these types of statistics.", type='str')

    parser.add_option('', '--single-sampler', dest='single_sampler', 
                      default=False, help='Use only a single sampler', action='store_true')
    parser.add_option('', '--no-rmsd', dest='no_rmsd', 
                      default=False, help='Refrain from trying to calculate the rmsd.', action='store_true')
    parser.add_option('', '--dists', dest='dists', default=None, 
                      help="Calculate the distance between pairs of nucleotides (i.e. 14,96:14,119)", 
                      type='str')
    parser.add_option('', '--save-iterative-cg-measures', dest='save_iterative_cg_measures', default=False, help='Save the coarse-grain measures every time the energy function is recalculated', action='store_true')
    parser.add_option('', '--jared-dir', dest='jared_dir', default=None, help='Use JAR3D to predict geometries for the interior loops', type='str')
    parser.add_option('', '--start-at-native', dest='start_at_native', default=False, action='store_true', help='Start at the native conformation')
    parser.add_option('', '--fix-loop', dest='fix_loop', default=None, help='Fix the correct coordinates of a particular loop to the correct ones')
    parser.add_option('', '--fix-all-loops', dest='fix_all_loops', default=False, action='store_true',  help='Fix the geometries of all loops in the structure')
    parser
    parser.add_option('', '--no-constraint', dest='no_constraint', default=False, action='store_true', help="Don't use a constraint energy")

    (options, args) = parser.parse_args()

    fud.pv('options.no_rmsd')


    if len(args) < 1:
        print "Usage: ./gibbs.py temp.comp"
        print "Or ./gibb.py temp.fa. If the extension of the argument file ends in .fa, then treat it as a fasta file."

        sys.exit(1)

    fud.pv('args')

    bgs = []


    for arg in args:
        if arg[-3:] == '.fa':
            bgs += bgs_from_fasta(arg)
        else:
            bgs += [ftmc.CoarseGrainRNA(arg)]

    if len(bgs) > 1:
        print >> sys.stderr, "WARNING: More than one structure entered... only the first will be bearbeitet"

    #bg.calc_bp_distances()

    energies_to_sample = []
    if options.cyl_intersect:
        energies_to_sample += [fbe.CombinedEnergy([], [fbe.CylinderIntersectionEnergy()])]
    if options.radius_of_gyration:
        sse = fbe.RadiusOfGyrationEnergy(energy_prefactor=options.energy_prefactor)
        sse.background = options.background
        energies_to_sample += [fbe.CombinedEnergy([], [sse])]

    if options.constant_energy:
        ce = fbe.ConstantEnergy()
        energies_to_sample += [fbe.CombinedEnergy([], [ce])]

    if options.random_energy:
        re = fbe.RandomEnergy()
        energies_to_sample += [fbe.CombinedEnergy([], [re])]

    if options.loop_energy:
        lle = fbe.ShortestLoopDistanceEnergy()
        energies_to_sample += [fbe.CombinedEnergy([], [lle])]

    if options.aminor_shortestloop:
        nonconstraint = []

        bg = bgs[0]
        nonconstraint += [fbe.ShortestLoopDistanceEnergy()]
        '''
        for hloop in bg.hloop_iterator():
            nonconstraint += [fbe.ShortestLoopDistancePerLoop(hloop)]
        '''

        nonconstraint += [fbe.AMinorEnergy(loop_type = 'h')]
        nonconstraint += [fbe.AMinorEnergy(loop_type = 'i')]

        energies_to_sample += [fbe.CombinedEnergy([], nonconstraint)]

    if options.specific_aminor:
        nonconstraint = [fbe.RadiusOfGyrationEnergy()]
        bg = bgs[0]

        # if we specify all, then we try and maximize the A-Minor interaction potential
        # for all internal and hairpin loops
        if len(options.specific_aminor.split(',')) == 1 and options.specific_aminor == 'all':
            for d in bg.defines:
                if d[0] == 'i' or d[0] == 'h':
                    if 'AA' in "".join(bg.get_define_seq_str(d)):
                        nonconstraint += [fbe.SpecificAMinorEnergy(loop_name=d, energy_prefactor=1)]
                        fud.pv('d')
        else:
            for sa in options.specific_aminor.split(','):
                nonconstraint += [fbe.SpecificAMinorEnergy(loop_name=sa, energy_prefactor=1)]

        for hloop in bg.hloop_iterator():
            if len(list(bg.define_residue_num_iterator(hloop))) > 4:
                fud.pv('hloop')
                nonconstraint += [fbe.ShortestLoopDistancePerLoop(hloop)]

        energies_to_sample += [fbe.CombinedEnergy([], nonconstraint)]

    if options.aminor_perloop:
        nonconstraint = []

        bg = bgs[0]
        for hloop in bg.hloop_iterator():
            nonconstraint += [fbe.ShortestLoopDistancePerLoop(hloop)]

        nonconstraint += [fbe.AMinorEnergy(loop_type = 'h')]
        nonconstraint += [fbe.AMinorEnergy(loop_type = 'i')]

        energies_to_sample += [fbe.CombinedEnergy([], nonconstraint)]

    if options.cylinder_perloop_radius_of_gyration:
        cie = fbe.CylinderIntersectionEnergy()
        #lle = fbe.ShortestLoopDistanceEnergy()
        rog = fbe.RadiusOfGyrationEnergy(energy_prefactor=options.energy_prefactor)
        nonconstraint = [rog, cie]

        bg = bgs[0]
        for hloop in bg.hloop_iterator():
            nonconstraint += [fbe.ShortestLoopDistancePerLoop(hloop)]

        rog.background = options.background
        energies_to_sample += [fbe.CombinedEnergy([], nonconstraint)]

    if options.cylinder_shortestloop_radius_of_gyration:
        cie = fbe.CylinderIntersectionEnergy()
        lle = fbe.ShortestLoopDistanceEnergy()
        rog = fbe.RadiusOfGyrationEnergy(energy_prefactor=options.energy_prefactor)
        rog.background = options.background
        energies_to_sample += [fbe.CombinedEnergy([], [lle, rog, cie])]

    if options.cylinder_loop_radius_of_gyration:
        cie = fbe.CylinderIntersectionEnergy()
        lle = fbe.LoopLoopEnergy()
        rog = fbe.RadiusOfGyrationEnergy(energy_prefactor=options.energy_prefactor)
        rog.background = options.background
        energies_to_sample += [fbe.CombinedEnergy([], [lle, rog, cie])]

    if options.cylinder_loop:
        lle = fbe.LoopLoopEnergy()
        cie = fbe.CylinderIntersectionEnergy()
        sse.background = options.background
        energies_to_sample += [fbe.CombinedEnergy([], [cie, lle])]

    if options.aminor:
        ame1 = fbe.AMinorEnergy(loop_type='h')
        ame2 = fbe.AMinorEnergy(loop_type='i')
        energies_to_sample += [fbe.CombinedEnergy([], [ame1, ame2])]
        #energies_to_sample += [fbe.CombinedEnergy([], [sse, ame1])]

    if options.aminor_radius_of_gyration:
        sse = fbe.RadiusOfGyrationEnergy(energy_prefactor=options.energy_prefactor)
        ame1 = fbe.AMinorEnergy(loop_type='h')
        ame2 = fbe.AMinorEnergy(loop_type='i')
        sse.background = options.background
        energies_to_sample += [fbe.CombinedEnergy([], [sse, ame1, ame2])]
        #energies_to_sample += [fbe.CombinedEnergy([], [sse, ame1])]

    if options.cylinder_radius_of_gyration:
        sse = fbe.RadiusOfGyrationEnergy(energy_prefactor=options.energy_prefactor)
        cie = fbe.CylinderIntersectionEnergy()
        sse.background = options.background
        energies_to_sample += [fbe.CombinedEnergy([], [sse, fbe.CylinderIntersectionEnergy()])]

    if options.loop_radius_of_gyration:
        sse = fbe.RadiusOfGyrationEnergy(energy_prefactor=options.energy_prefactor)
        sse.background = options.background
        energies_to_sample += [fbe.CombinedEnergy([], [sse, fbe.LoopLoopEnergy()])]

    if options.distance_energy:
        energies_to_sample += [fbe.DistanceEnergy(bg.get_long_range_constraints())]


    if len(energies_to_sample) == 0 or options.aminor_perloop_radius_of_gyration:
        rog = fbe.RadiusOfGyrationEnergy(energy_prefactor=options.energy_prefactor)
        nonconstraint = [rog]

        bg = bgs[0]
        for hloop in bg.hloop_iterator():
            nonconstraint += [fbe.ShortestLoopDistancePerLoop(hloop)]
        nonconstraint += [fbe.AMinorEnergy(loop_type = 'h')]
        nonconstraint += [fbe.AMinorEnergy(loop_type = 'i')]
        nonconstraint += [fbe.StemVirtualResClashEnergy()]

        rog.background = options.background
        energies_to_sample += [fbe.CombinedEnergy([], nonconstraint)]

    if options.clamp is not None:
        pairs = options.clamp.split(':')
        bg = bgs[0]


        for p in pairs:
            r1,r2 = p.split(',')

            try:
                # initially we assume the clamp target are residue numbers
                r1 = int(r1)
                r2 = int(r2)

                e1 = bg.get_node_from_residue_num(int(r1))
                e2 = bg.get_node_from_residue_num(int(r2))
            except ValueError:
                # ... or they are element names
                e1 = r1
                e2 = r2

            if e1 not in bg.defines.keys() or e2 not in bg.defines.keys():
                print >>sys.stderr, "ERROR: Invalid values for clamping"
                sys.exit(1)

            if e1 == e2:
                print >>sys.stderr, "Can't clamp identical elements"

            print >>sys.stderr, "clamping {0}, {1}".format(e1,e2)
            # the first energy to sample should be a CombinedEnergy
            energies_to_sample[0].energies += [fbe.DistanceExponentialEnergy(e1,e2,15.,1.)]

    for bg in bgs:
        options.bg_filename = args[0]
        fud.pv('energies_to_sample')

        if len(list(bg.stem_iterator())) == 0:
            print >> sys.stderr, "Cannot simulate an open chain, the structure needs to have at least one stem"
            sys.exit(1)

        predict(bg, energies_to_sample, options)

if __name__ == '__main__':
    seed_num = random.randint(0, sys.maxint)
    try:
        random.seed(seed_num)
        main()
    except:
        #pdb.post_mortem()
        print >>sys.stderr, "random seed:", seed_num
        raise

