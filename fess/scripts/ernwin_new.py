#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip) #future package
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)

import argparse, sys, warnings, copy, os, random, math
import os.path
import subprocess as spr
import contextlib
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.stats as ftms
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.utilities.debug as fud
from fess.builder import energy as fbe
from fess.builder import models as fbm
from fess.builder import sampling as fbs
from fess.builder import config
from fess.builder import samplingStatisticsNew2 as sstats
from fess import data_file
import scipy.ndimage
#Magic numbers
DEFAULT_ENERGY_PREFACTOR=30


@contextlib.contextmanager
def open_for_out(filename=None):
    if filename and filename != '-':
        fh = open(filename, 'w')
    else:
        fh = sys.stdout
    try:
        yield fh
    finally:
        if fh is not sys.stdout:
            fh.close()

def get_parser():
    """
    Here all commandline & help-messages arguments are defined.

    :returns: an instance of argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
    #Argument(s)
    parser.add_argument('rna', nargs=1, help='A *.fa or *.cg/*.coord file holding the RNA to analyze.')
    #Options
    #Modify general behavior
    parser.add_argument('-i', '--iterations', action='store', default=10000, help='Number of structures to generate', type=int)
    parser.add_argument('--start-from-scratch', default=False, action='store_true', 
                        help="Do not attempt to start at the input conformation.\n"
                             "(Automatically True for fasta files.)")
    parser.add_argument('--eval-energy', default=False, action='store_true', 
                        help='Evaluate the energy of the input structure and\n'
                             'exit without sampling.')
    parser.add_argument('--seed', action='store', help="Seed for the random number generator.",
                        type=int)
    parser.add_argument('--exhaustive', type=str, help="A STRING describing the coarse grained element (e.g. i2)\n"
                                                       "Instead of the MCMC sampling algorithm, use a sampler \n"
                                                       "that tries every stat for the given coarse grain element.")
    parser.add_argument('--new-ml', action="store_true", help="Experimental")
    #Controll output
    parser.add_argument('--save-n-best', default=3, 
                        help='Save the best (lowest energy) n structures.', type=int)
    parser.add_argument('--save-min-rmsd', default=3, 
                        help='Save the best (lowest rmsd) n structures.', type=int)
    parser.add_argument('--step-save', default=0, help="Save the structure at every n'th step.",
                         type=int)
    parser.add_argument('--dump-energies', default=False, action='store_true', 
                        help='Dump the measures used for energy calculation to file') #UNUSED OPTION. REMOVE
    parser.add_argument('--no-rmsd', default=False, 
                        help='Refrain from trying to calculate the rmsd.', action='store_true')
    parser.add_argument('--rmsd-to', action='store', type=str, 
                        help="A *.cg/ *.coord or *.pdb file.\n"
                             "Calculate the RMSD and MCC relative to the structure in this file,\n"
                             "not to the structure used as starting point for sampling.")
    #Controll output files
    parser.add_argument('--output-file', action='store', type=str, default="out.log",
                        help="Filename for output (log). \n"
                             "This file will be created inside the --output-base-dir.\n"
                             "Default: out.log")
    parser.add_argument('--output-base-dir', action='store', type=str, default="", 
                        help="Base dir for the output. \n"
                              "In this directory, a subfolder with the name\n"
                              "from the fasta-file will be generated")
    parser.add_argument('-o', '--output-dir-suffix', action='store', type=str, 
                        default="", help="Suffix attached to the name from the fasta-file, \n"
                                         "used as name for the subfolder with all structures.")
    #Controll Stats for sampling
    parser.add_argument('--stats-file', type=str, default=data_file("stats/all.stats"),
                        help= "A filename.\n"
                              "A file containing all the stats to sample from\n"
                              " for all coarse grained elements")
    parser.add_argument('--clustered-angle-stats', type=str, action="store",
                        help= "A filename.\n"
                              "If given, use this instead of --stats-file for the\n"
                              " angle stats (interior and multi loops).\n" 
                              "A clustered stats file can be created with scrips/cluster_stats.py\n"
                              "This is used to sample equally from all CLUSTERS of angle stats.\n"
                              "And is used to compensate for unequally populated clusters.\n"
                              "The statistical bias introduced by this file can be compensated\n"
                              " for with the SBC (Stats bias compensation energy) [TODO]")
    parser.add_argument('--jar3d-dir', type=str, help="The base dir of your JAR3D (Motiv atlas) installation.\n"
                                                      "It should contain the 'JAR3D' subdirectory and \n"
                                                      "the file 'scripts/annotate_structure.py'\n"
                                                      "JAR3D is available at 'https://github.com/BGSU-RNA/JAR3D'\n"
                                                      "or 'http://rna.bgsu.edu/jar3d'")
    #Choose energy function(s)
    parser.add_argument('-c', '--constraint-energy', default="D", action='store', type=str, 
                                    help="The type of constraint energy to use. \n"
                                         "D=Default    clash- and junction closure energy\n"
                                         "B=Both       same as 'D'\n"
                                         "J=junction   only junction closure energy\n"
                                         "C=clash      only stem clash energy")
    parser.add_argument('-e', '--energy', default="D", action='store', type=str, 
                        help= "The type of non-constraint energy to use. D=Default\n"
                              "Specify a ','-separated list of energy contributions.\n"
                              "Each contribution has the format: [PRE]TYP[ADJ].\n"
                              "PRE: optional energy prefactor\n"
                              "ADJ: optional adjustment of target distribution\n"
                              "     (float), default=1.0\n"
                              "     use START_END or START_STEP_END to modify the\n"
                              "     adjustment during sampling. E.g. 1.0_0.1_1.4\n"
                              "     For each adjustment, equally many sampling steps \n"
                              "     are used.\n"
                              "     If step is not given, 1.0 or 0.1 is used, \n"
                              "     depending on the difference between START and END.\n"
                              "TYP: one of the following\n"
                              "       ROG:  Radius of gyration energy\n"
                              "       NDR:  Normal Distributed Radius of Gyration energy.\n"
                              "             Use a normal distribution for the target ROG\n"
                              "             With 0.77*ADJ as mean and 0.23*ADJ stddev\n."
                              "             [0.77 is a rough estimate for the relation between\n"
                              "             perimeter and radius of gyration]\n"
                              "       SLD:  shortest loop distance per loop\n"
                              "       AME:  A-Minor energy\n"
                              "       PRO:  Match Projection distances. \n"
                              "             Requires the --projected-dist option\n"
                              "       HDE:  Hausdorff distance based Energy \n"
                              "             Requires the --ref-img and --scale  options.\n"
                              "       DEF:  add all default energies to the \n"
                              "             combined energy\n"
                              "       CHE:  Cheating Energy. Tries to minimize the RMSD\n"
                              "       CLA:  Clamp elements (not nucleotides) together \n"
                              "             (at 15 Angstrom) with an exponential energy. \n"
                              "             The prefactor is used as scale (default 1).\n"
                              "             Requires the --clamp option\n"
                              "Example: ROG10,SLD,AME")
    parser.add_argument('--track-energies', action='store', type=str, default="",
                        help= "A ':' seperated list of combined energies.\n"
                              "Each energy is given in the format specified\n"
                              "for the --energy option.\n"
                              "These energies are not used for sampling, \n"
                              "but only calculated after each step.")
    parser.add_argument('--projected-dist', action='store', type=str, default="",
                        help= "A ':' seperated list of tripels: \n"
                              "cgelement, cgelement, dist\n"
                              "Where dist is the projected distance in the image \n"
                              "in Angstrom.\n"
                              "Example: 's1,h3,10:s2,h3,12'")
    parser.add_argument('--ref-img', action='store', type=str, default="",
                        help= "A black and white square image (e.g. in png format)\n"
                              "as a reference projection for the Hausdorff Energy.\n"
                              "White is the RNA, black is the background.\n"
                              "Requires the Python Imaging Library (PIL) or Pillow.")
    parser.add_argument('--scale', action='store', type=int,
                        help= "Used for the Hausdorff Energy.\n"
                              "The length (in Angstrom) of each side \n"
                              "of the image is")
    parser.add_argument('--clamp', action='store', type=str,
                        help= "Used for the CLA energy.\n"
                              "A list `p1,p2:p3,p4:...` where p1 and p2 are clamped together\n"
                              " and p3+p4 are clamped together. The pi are either emelents\n"
                              " ('s1','i1',...) or integers (positions in the sequence).\n")
    return parser

def getSLDenergies(cg, prefactor=DEFAULT_ENERGY_PREFACTOR):
    """
    Get the shortest loopdistance per loop energy for each hloop.

    :param cg: The coarse grained RNA
    :returns: A list of energies
    """
    energies=[]
    for hloop in cg.hloop_iterator():
        energies+= [fbe.ShortestLoopDistancePerLoop(hloop, prefactor)]
    return energies
def getAMinorEnergies(cg, pre=DEFAULT_ENERGY_PREFACTOR, adj=1.0):
    """
    Get the A-minor energies for h- and i-loops.
    
    :param pre: Energy prefactor
    :param adj: Adjustment

    :returns: A list of energies
    """
    ame1 = fbe.AMinorEnergy(loop_type = 'h', energy_prefactor=pre, adjustment=adj)
    ame2 = fbe.AMinorEnergy(loop_type = 'i', energy_prefactor=pre, adjustment=adj)
    energies = []
    if ame1.get_num_loops(cg)>0:
        energies.append(ame1)
    if ame2.get_num_loops(cg)>0:
        energies.append(ame2)
    return energies

def getDefaultEnergies(cg):
    """
    Get the default energies.

    :param cg: The coarse grained RNA
    :returns: A list of energies
    """
    energies=[]
    #Radius of gyration
    energies.append(fbe.RadiusOfGyrationEnergy(energy_prefactor=DEFAULT_ENERGY_PREFACTOR,
                                               adjustment=1.0))
    #Shortest loop distance per loop
    sld = getSLDenergies(cg)
    if sld:
        energies+=[ fbe.CombinedEnergy([],  sld , normalize=True) ]
    #A-minor energies
    energies += getAMinorEnergies(cg)
    return energies

def getPROenergy(proj_dist_str, prefactor):
    contributions=proj_dist_str.split(":")
    distances={}
    for contrib in contributions:
        f=contrib.split(",")
        if len(f)!=3:
            raise ValueError("Could not parse projected-dist string: '{}'".format(contrib))
        distances[(f[0],f[1])]=float(f[2])
    if not distances:
        raise ValueError("The --projected-dist option is required if the PRO energy is used.".format(contrib))
    return fbe.ProjectionMatchEnergy(distances, prefactor)

def getHDEenergy(hde_image, scale, pre):
    img=scipy.ndimage.imread(hde_image)
    return fbe.HausdorffEnergy(img, scale, pre)

def parseCombinedEnergyString(stri, cg, iterations, proj_dist, scale, hde_image, reference_cg):
    """
    Parses an energy string, as used for the --energy commandline option

    :param stri: The energy string
    :param cg: The coarse-grained RNA
    :iterations: The number of iterations. This is used for adjustments that change during sampling
    :returs: A combined energy
    """
    contributions=stri.split(",")
    energies=[]
    for contrib in contributions:
        if "DEF" in contrib:
            pre,_,adj=contrib.partition("DEF")
            if pre!="" or adj!="":
                warnings.warn("Prefactor '{}' and adjustment '{}' are ignored "
                              "for default energy!".format(pre, adj))
            energies+=getDefaultEnergies(cg)
        elif "ROG" in contrib:
            pre, adj=parseEnergyContributionString(contrib, "ROG")
            e = fbe.RadiusOfGyrationEnergy( energy_prefactor=pre, adjustment=adj[0])
            if adj[1]:
                e.set_dynamic_adjustment(adj[1],math.ceil(iterations/adj[2]))
            energies.append(e)
        elif "NDR" in contrib:
            pre, adj=parseEnergyContributionString(contrib, "NDR")
            e = fbe.NormalDistributedRogEnergy( energy_prefactor=pre, adjustment=adj[0])
            if adj[1]:
                e.set_dynamic_adjustment(adj[1],math.ceil(iterations/adj[2]))
            energies.append(e)
        elif "AME" in contrib:
            pre, adj=parseEnergyContributionString(contrib, "AME")
            es=getAMinorEnergies(cg, pre, adj[0])
            if adj[1]:
                for e in es:
                    e.set_dynamic_adjustment(adj[1],math.ceil(iterations/adj[2]))
            energies+=es
        elif "SLD" in contrib:
            pre,_, adj=contrib.partition("SLD")
            if adj!="":
                warnings.warn("Adjustment '{}' is ignored for "
                              "ShortestLoopdistancePerLoop energy!".format(adj))            
            if pre=="": pre=DEFAULT_ENERGY_PREFACTOR
            else: pre=int(pre)
            energies+=[ fbe.CombinedEnergy([],  getSLDenergies(cg, pre), normalize=True ) ]
        elif "PRO" in contrib:
            pre,_, adj=contrib.partition("PRO")
            if adj!="":
                warnings.warn("Adjustment '{}' is ignored for "
                              "ProjectionMatchEnergy energy!".format(adj))
            if pre:
                pre=int(pre)
            else:
                pre=1
            energies.append(getPROenergy(proj_dist, pre))
        elif "HDE" in contrib:
            pre,_, adj=contrib.partition("HDE")
            if adj!="":
                warnings.warn("Adjustment '{}' is ignored for Hausdorff energy!".format(adj))
            if pre:
                pre=int(pre)
            else:
                pre=DEFAULT_ENERGY_PREFACTOR
            energies.append(getHDEenergy(hde_image, scale, pre))
        elif "CHE" in contrib:
            pre,_, adj=contrib.partition("CHE")
            if pre!="" or adj!="":
                warnings.warn("Prefactor '{}' and adjustment '{}' are ignored "
                              "for cheating energy!".format(pre, adj))
            energies.append(fbe.CheatingEnergy(reference_cg))
        elif "CLA" in contrib:
            pre,_, adj=contrib.partition("CLA")
            if adj!="":
                warnings.warn("Adjustment '{}' is ignored for Clamp Energy!".format(adj))
            if pre:
                pre=int(pre)
            else:
                pre=1
            pairs = args.clamp.split(':')
            clamp=[]
            for pair in pairs:
                r1,r2=pair.split(",")
                try: # initially we assume the clamp target are residue numbers
                    r1=int(r1)
                except ValueError: #Or they are element names
                    e1=r1
                else:
                    e1=cg.get_node_from_residue_num(r1)
                try: # initially we assume the clamp target are residue numbers
                    r2=int(r2)
                except ValueError: #Or they are element names
                    e2=r2
                else:
                    e2=cg.get_node_from_residue_num(r2)
                if e1 not in cg.defines or e2 not in cg.defines:
                    print("ERROR: Invalid Clamp values '{}'-'{}' "
                          "(elements '{}'-'{}').".format(r1,r2,e1,e2), file=sys.stderr)
                    sys.exit(1)
                if e1==e2:
                    warnings.warn("Cannot clamp identical elements "
                                  " {} and {} ({}=={})".format(r1,r2,e1,e2))
                else:
                    clamp+=[fbe.DistanceExponentialEnergy(e1,e2,15.,pre)]
            if clamp:
                energies.append(fbe.CombinedEnergy([],clamp))
        else:
            print("ERROR: Cannot parse energy contribution: '{}'".format(contrib), file=sys.stderr)
            sys.exit(1)
    return fbe.CombinedEnergy([], energies)
def parseEnergyContributionString(contrib, sep):
    """
    Helper function for parsing of commandline.

    Splits the contrib string at the occurence of sep. 
    If an empty string is before and/por after the seperator, the default value is returned.

    :param contrib: A string, e.g. "30ROG1.0"
    :returns: a 2-tuple of floats: Prefactor, Adjustment
    """
    pre,_,adj=contrib.partition(sep)
    if not pre:
        pre=DEFAULT_ENERGY_PREFACTOR
    if "_" in adj:
        a=adj.split("_")
        if len(a)==2:
            adj=float(a[0])
            end=float(a[1])
            if abs(adj-end)<1.2:
                step=0.1
            else:
                step=1
            if end<adj:
                step=step*-1
            num_steps=math.ceil((end-adj)/step)+1
            assert num_steps>1, numSteps
            return float(pre), (adj, step, num_steps)
        elif len(a)==3:
            adj=float(a[0])
            step=float(a[1])
            end=float(a[2])
            num_steps=math.ceil((end-adj)/step)+1
            if num_steps<2:
                raise ValueError("Could not parse adjustment in option '{}': "
                                  "Expected START_STEP_STOP, found {}, {}, {} which would "
                                  "lead to a total of {} steps".format(contrib, adj, step, 
                                                                       end, num_steps))
            return float(pre), (adj, step, num_steps)
        else:
            raise ValueError("Could not parse adjustment in option: '{}'".format(contrib))
        
    if not adj:
        adj=1.0
    return float(pre), (float(adj),0,0)



def setup_deterministic(args):
    """
    The part of the setup procedure that does not use any call to the random number generator.

    :param args: An argparse.ArgumentParser object holding the parsed arguments.
    """
    #Load the RNA from file
    rnafile, = args.rna #Tuple unpacking
    if rnafile[-3:] == '.fa':
        cg = ftmc.CoarseGrainRNA()
        try:
            with open(rnafile) as fastafile:
                cg.from_fasta(fastafile.read())
        except IOError as e:
            print("ERROR: Could not open file '{}' for reading:".format(rnafile),e, file=sys.stderr)
            sys.exit(1)
    else:
        cg = ftmc.CoarseGrainRNA(rnafile)

    #Output file and directory        
    ofilename=None
    if not args.eval_energy:
        if args.output_base_dir and not os.path.exists(args.output_base_dir):
            os.makedirs(args.output_base_dir)
            print ("INFO: Directory {} created.".format(args.output_base_dir), file=sys.stderr)
        subdir=cg.name+args.output_dir_suffix
        config.Configuration.sampling_output_dir = os.path.join(args.output_base_dir, subdir)
        if not os.path.exists(config.Configuration.sampling_output_dir):
            os.makedirs(config.Configuration.sampling_output_dir)
            print ("INFO: Directory {} created. This folder will be used for all output "
                   "files.".format(config.Configuration.sampling_output_dir), file=sys.stderr)
        if args.output_file:
            ofilename=os.path.join(config.Configuration.sampling_output_dir, args.output_file)


    #Initialize the spatial model
    if args.clustered_angle_stats and args.jar3d_dir:
        print("ERROR: --clustered-angle-stats and --jar3d-dir are mutually exclusive!", file=sys.stderr)
        sys.exit(1)
    if args.clustered_angle_stats:
        sm=fbm.SpatialModel(cg, ftms.get_conformation_stats(args.stats_file, args.clustered_angle_stats))
    elif args.jar3d_dir:
        jared_script = op.join(options.jared_dir, 'scripts/annotate_structure.py')
        jared_data   = op.join(options.jared_dir, 'JAR3D')
        jared_out    = op.join(config.Configuration.sampling_output_dir, filtered_stats)
        
        cmd = ['python', jared_script, options.bg_filename, '-o', jared_data,
               '-m', '-e', '-d', jared_data]

        p = spr.Popen(cmd, stdout=spr.PIPE)
        out, err = p.communicate()

        with open(jared_out, 'w') as filtered_out:
            filtered_out.write(out)

        filtered_stats = ftms.FilteredConformationStats(stats_file=args.stats_file,
                                                        filter_filename=jared_out)
        ftms.set_conformation_stats(filtered_stats)
        sm=fbm.SpatialModel(cg, ftms.get_conformation_stats())
    else:
        sm=fbm.SpatialModel(cg, ftms.get_conformation_stats(args.stats_file))
    #Load the reference sm (if given)
    if args.rmsd_to:
        if args.rmsd_to.endswith(".pdb"):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pdb = list(bp.PDBParser().get_structure('reference', args.rmsd_to).get_chains())[0]
                original_sm  = fbm.SpatialModel(ftmc.from_pdb(pdb))
        else:
            original_sm=fbm.SpatialModel(ftmc.CoarseGrainRNA(args.rmsd_to))
    else:
        original_sm=fbm.SpatialModel(copy.deepcopy(sm.bg))
        


    #Initialize the requested energies
    
    if args.energy=="D":
        energy=fbe.CombinedEnergy([],getDefaultEnergies(cg))     
    else:
        energy=parseCombinedEnergyString(args.energy, cg, args.iterations,
                                         args.projected_dist,args.scale,
                                         args.ref_img, original_sm.bg)

    #Initialize energies to track
    energies_to_track=[]
    for track_energy_string in args.track_energies.split(":"):
        if track_energy_string:
            if track_energy_string=="D":
                energies_to_track.append(fbe.CombinedEnergy([],getDefaultEnergies(cg)))
            else:
                energies_to_track.append(parseCombinedEnergyString(track_energy_string, cg,
                                         args.iterations, args.projected_dist, args.scale,
                                         args.ref_img, original_sm.bg))


    #Initialize the Constraint energies
    junction_energy=None
    clash_energy=None
    if args.constraint_energy in ["D","B","J"]:
        sm.junction_constraint_energy=fbe.CombinedEnergy([fbe.RoughJunctionClosureEnergy()])
    if args.constraint_energy in ["D","B","C"]:
        sm.constraint_energy=fbe.CombinedEnergy([fbe.StemVirtualResClashEnergy()])
    return sm, original_sm, ofilename, energy, energies_to_track

def setup_stat(out_file, sm, args, energies_to_track, original_sm):
    """
    Setup the stat object used for logging/ output.

    :param out_file: an opened file handle for writing
    :param sm: The spatial model. Note: A deepcopy of this object will be generated as a 
               reference structure. This is unused if args.rmds_to is set.
    :param args: An argparse.ArgumentParser object holding the parsed arguments.
    """


    #stat = fbs.SamplingStatistics(original_sm, plter , None, silent=False,
    #                                  output_file=out_file, 
    #                                  save_n_best = args.save_n_best, 
    #                                  dists = [], 
    #                                  save_iterative_cg_measures=args.save_iterative_cg_measures, 
    #                                  no_rmsd = args.no_rmsd)
    #stat.step_save = args.step_save
    options={}
    if args.no_rmsd:
        options["rmsd"] = False
        options["acc"]  = False
    if not args.start_from_scratch and not args.rmsd_to:
        options["extreme_rmsd"] = "max" #We start at 0 RMSD. Saving the min RMSD is useless.
    options[ "step_save" ] = args.step_save 
    options[ "save_n_best" ] = args.save_n_best 
    options[ "save_min_rmsd" ] = args.save_min_rmsd
    options[ "measure" ]=[]
    for energy in energies_to_track:
        if  (isinstance(energy, fbe.ProjectionMatchEnergy) or
            isinstance(energy, fbe.HausdorffEnergy)):
                options[ "measure" ].append(energy)
    stat = sstats.SamplingStatistics(original_sm, energy_functions = energies_to_track,
                                     output_file=out_file, options=options)
    return stat


# Parser is available even if __name__!="__main__", to allow for 
# documentation with sphinxcontrib.autoprogram
parser = get_parser()
if __name__=="__main__":
    args = parser.parse_args()

    #Setup that does not use the random number generator.
    randstate=random.getstate()#Just for verification purposes
    sm, original_sm, ofilename, energy, energies_to_track = setup_deterministic(args)
    assert randstate==random.getstate()#Just for verification purposes
    fud.pv("energies_to_track")
    #Eval-energy mode
    if args.eval_energy:
        sm.bg.add_all_virtual_residues()
        fud.pv('energy.eval_energy(sm, verbose=True, background=False)')
        if sm.constraint_energy:
            fud.pv('sm.constraint_energy.eval_energy(sm, verbose=True, background=False)')
        if sm.junction_constraint_energy:
            fud.pv('sm.junction_constraint_energy.eval_energy(sm, verbose=True, background=False)')
        for track_energy in energies_to_track:
            fud.pv('track_energy.eval_energy(sm, verbose=True, background=False)')
        sys.exit(0) 
  
    #Set-up the random Number generator.
    #Until here, no call to random should be made.
    if args.seed:
        seed_num=args.seed
    else:
        seed_num = random.randint(0, sys.maxint)            
    random.seed(seed_num)
    #Main function, dependent on random.seed        
    plter=None #fbs.StatisticsPlotter()
    with open_for_out(ofilename) as out_file:
        if isinstance(energy, fbe.CombinedEnergy):
            energies_to_track+=energy.uncalibrated_energies
        elif isinstance(energy, fbe.CoarseGrainEnergy):
            energies_to_track+=[energy]
        stat=setup_stat(out_file, sm, args, energies_to_track, original_sm)
        try:
            print ("# Random Seed: {}".format(seed_num), file=out_file)
            print ("# Command: `{}`".format(" ".join(sys.argv)), file=out_file)
            if args.exhaustive:
                sampler = fbs.ExhaustiveExplorer(sm, energy, stat, args.exhaustive, args.start_from_scratch)
            elif args.new_ml:
                sampler = fbs.ImprovedMultiloopMCMC(sm, energy, stat, 
                                          start_from_scratch=args.start_from_scratch,
                                          dump_measures=args.dump_energies)
            else:
                sampler = fbs.MCMCSampler(sm, energy, stat, 
                                          start_from_scratch=args.start_from_scratch,
                                          dump_measures=args.dump_energies)
            for i in range(args.iterations):
                sampler.step()            
            #plter.finish()
        finally: #Clean-up 
            print("INFO: Random seed was {}".format(seed_num), file=sys.stderr)
