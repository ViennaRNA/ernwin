#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip)
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)

import argparse, sys, warnings, copy, os, random
import os.path
import contextlib
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.utilities.debug as fud
from fess.builder import energy as fbe
from fess.builder import models as fbm
from fess.builder import sampling as fbs
from fess.builder import config
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
    parser.add_argument('--start-from-scratch', default=False, action='store_true', help="Do not attempt to start at the input conformation.\n"
                                                                                         "(Automatically True for fasta files.)")
    parser.add_argument('--eval-energy', default=False, action='store_true', help='Evaluate the energy of the input structure and\n'
                                                                                  'exit without sampling.')
    parser.add_argument('--seed', action='store', help="Seed for the random number generator.", type=int)
    #Controll output
    parser.add_argument('--save-n-best', default=3, help='Save the best n structures.', type=int)
    parser.add_argument('--step-save', default=0, help="Save the structure at every n'th step.", type=int)
    parser.add_argument('--save-iterative-cg-measures', default=False, help="Save the coarse-grain measures every time\n"
                               " the energy function is recalculated", action='store_true')
    parser.add_argument('--dump-energies', default=False, action='store_true', help='Dump the energies to file')
    parser.add_argument('--no-rmsd', default=False, help='Refrain from trying to calculate the rmsd.', action='store_true')
    #Controll output files
    parser.add_argument('--output-file', action='store', type=str, help="Filename for output (log). \n"
                                                    "This file will be created inside the --output-base-dir.\n"
                                                    "Default: standard out")
    parser.add_argument('--output-base-dir', action='store', type=str, default="", help="Base dir for the output. \n"
                                                   "In this directory, a subfolder with the name\n"
                                                   "from the fasta-file will be generated")
    parser.add_argument('-o', '--output-dir-suffix', action='store', type=str, default="", help="Suffix attached to the name from the fasta-file, \n"
                                                                                                "used as name for the subfolder with all structures.")
    #Choose energy function(s)
    parser.add_argument('-c', '--constraint-energy', default="D", action='store', type=str, help="The type of constraint energy to use. \n"
                                                   "D=Default    clash- and junction closure energy\n"
                                                   "B=Both       same ad 'D'\n"
                                                   "J=junction   only junction closure energy\n"
                                                   "C=clash      only stem clash energy")
    parser.add_argument('-e', '--energy', default="D", action='store', type=str, 
                        help= "The type of non-constraint energy to use. D=Default\n"
                              "Specify a ','-separated list of energy contributions.\n"
                              "Each contribution has the format: [PRE]TYP[ADJ].\n"
                              "PRE: optional energy prefactor\n"
                              "ADJ: optional adjustment of target distribution\n"
                              "     (float), default=1.0\n"
                              "TYP: one of the following\n"
                              "       ROG:  Radius of gyration energy\n"
                              "       SLD:  shortest loop distance per loop\n"
                              "       AME:  A-Minor energy\n"
                              "       DEF:  add all default energies to the \n"
                              "             combined energy\n"
                              "Example: ROG10,SLD,AME")
    parser.add_argument('--track-energies', action='store', type=str, default="",
                        help= "A ':' seperated list of combined energies.\n"
                              "Each energy is given in the format specified\n"
                              "for the --energy option.\n"
                              "These energies are not used for sampling, \n"
                              "but only calculated after each step.")
    return parser

def getSLDenergies(cg):
    """
    Get the shortest loopdistance per loop energy for each hloop.

    :param cg: The coarse grained RNA
    :returns: A list of energies
    """
    energies=[]
    for hloop in cg.hloop_iterator():
        energies+= [fbe.ShortestLoopDistancePerLoop(hloop)]
    return energies
def getAMinorEnergies(pre=DEFAULT_ENERGY_PREFACTOR, adj=1.0):
    """
    Get the A-minor energies for h- and i-loops.
    
    :param pre: Energy prefactor
    :param adj: Adjustment

    :returns: A list of energies
    """
    return [fbe.AMinorEnergy(loop_type = 'h', energy_prefactor=pre, adjustment=adj),
            fbe.AMinorEnergy(loop_type = 'i', energy_prefactor=pre, adjustment=adj)]

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
    energies+=getSLDenergies(cg)
    #A-minor energies
    energies += getAMinorEnergies()
    return energies

def parseCombinedEnergyString(stri, cg):
    """
    Parses an energy string, as used for the --energy commandline option

    :param stri: The energy string
    :param cg: The coarse-grained RNA
    :returs: A combined energy
    """
    contributions=stri.split(",")
    energies=[]
    for contrib in contributions:
        if "DEF" in contrib:
            pre,_,adj=contrib.partition("DEF")
            if pre!="" or adj!="":
                warnings.warn("Prefactor '{}' and adjustment '{}' are ignored for default energy!".format(pre, adj))
            energies+=getDefaultEnergies(cg)
        elif "ROG" in contrib:
            pre, adj=parseEnergyContributionString(contrib, "ROG")
            energies.append(fbe.RadiusOfGyrationEnergy( energy_prefactor=pre, adjustment=adj))
        elif "AME" in contrib:
            pre, adj=parseEnergyContributionString(contrib, "AME")
            energies+=getAMinorEnergies(pre, adj)
        elif "SLD" in contrib:
            pre,_, adj=contrib.partition("SLD")
            if pre!="" or adj!="":
                warnings.warn("Prefactor '{}' and adjustment '{}' are ignored for ShortestLoopdistancePerLoop energy!".format(pre, adj))
            energies+=getSLDenergies(cg)
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
    if not adj:
        adj=1.0
    return float(pre), float(adj)



def setup_deterministic(args):
    """
    The part of the setup procedure that does not use any call to the random number generator.

    @TESTS: Verify that subsequent calls to this function with the same input lead to the same output. (TODO)

    :param args: An argparse.ArgumentParser object holding the parsed arguments.
    """
    #Load the RNA from file
    rnafile, = args.rna #Tuple unpacking
    if rnafile[-3:] == '.fa':
        cg = ftmc.CoarseGrainRNA()
        try:
            with open(rnafile) as fastafile:
                cg.from_fasta(rnafile)
        except IOError as e:
            print("ERROR: Could not open file '{}' for reading:".format(rnafile),e, file=sys.stderr)
            sys.exit(1)
    else:
        cg = ftmc.CoarseGrainRNA(rnafile)

    #Output file and directory
    if args.output_base_dir and not os.path.exists(args.output_base_dir):
        os.makedirs(args.output_base_dir)
        print ("INFO: Directory {} created.".format(args.output_base_dir), file=sys.stderr)
    subdir=cg.name+args.output_dir_suffix
    config.Configuration.sampling_output_dir = os.path.join(args.output_base_dir, subdir)
    if not os.path.exists(config.Configuration.sampling_output_dir):
        os.makedirs(cbc.Configuration.sampling_output_dir)
        print ("INFO: Directory {} created. This folder will be used for all output files.".format(cbc.Configuration.sampling_output_dir), file=sys.stderr)
    ofilename=None
    if args.output_file:
        ofilename=os.path.join(config.Configuration.sampling_output_dir, args.output_file)

    #Initialize the requested energies
    if args.energy=="D":
        energy=fbe.CombinedEnergy([],getDefaultEnergies(cg))     
    else:
        energy=parseCombinedEnergyString(args.energy, cg)
    #Initialize energies to track
    energies_to_track=[]
    for track_energy_string in args.track_energies.split(":"):
        if track_energy_string:
            if track_energy_string=="D":
                energies_to_track.append(fbe.CombinedEnergy([],getDefaultEnergies(cg)))
            else:
                energies_to_track.append(parseCombinedEnergyString(track_energy_string, cg))
        #Initialize the spatial model
    sm=fbm.SpatialModel(cg)

    #Initialize the Constraint energies
    junction_energy=None
    clash_energy=None
    if args.constraint_energy in ["D","B","J"]:
        sm.junction_constraint_energy=fbe.CombinedEnergy([fbe.RoughJunctionClosureEnergy()])
    if args.constraint_energy in ["D","B","C"]:
        sm.constraint_energy=fbe.CombinedEnergy([fbe.StemVirtualResClashEnergy()])
    



    return sm, ofilename, energy, energies_to_track

def setup_stat(out_file, sm, args):
    """
    Setup the stat object used for loging/ output.

    :param out_file: an opened file handle for writing
    :param sm: The spatial model. Note: A deepcopy of this object will be generated as a reference structure.
    :param args: An argparse.ArgumentParser object holding the parsed arguments.
    """
    original_sm=copy.deepcopy(sm)
    stat = fbs.SamplingStatistics(original_sm, None, None, silent=False,
                                      output_file=out_file, 
                                      save_n_best = args.save_n_best, 
                                      dists = [], 
                                      save_iterative_cg_measures=args.save_iterative_cg_measures, 
                                      no_rmsd = args.no_rmsd)
    stat.step_save = args.step_save
    return stat


if __name__=="__main__":
    parser = get_parser()
    args = parser.parse_args()

    #Setup that does not use the random number generator.
    randstate=random.getstate()#Just for verification purposes
    sm, ofilename, energy, energies_to_track = setup_deterministic(args)
    assert randstate==random.getstate()#Just for verification purposes
    fud.pv("energies_to_track")
    #Eval-energy mode
    if args.eval_energy:
        for s in sm.bg.stem_iterator():
            ftug.add_virtual_residues(sm.bg, s)
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
    with open_for_out(ofilename) as out_file:
        stat=setup_stat(out_file, sm,args)            
        try:
            print ("# Random Seed: {}".format(seed_num), file=out_file)
            sampler = fbs.MCMCSampler(sm, energy, stat, args.no_rmsd, energies_to_track=energies_to_track, start_from_scratch=args.start_from_scratch)
            sampler.dump_measures = args.dump_energies
            for i in range(args.iterations):
                sampler.step()
        finally: #Clean-up 
            print("INFO: Random seed was {}".format(seed_num), file=sys.stderr)
