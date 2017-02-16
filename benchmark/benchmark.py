#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip) #future package
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)
import matplotlib as mpl
mpl.use('Agg')
import logging    
logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)
import argparse
import glob
import os.path as op
import re
import numpy as np
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.ensemble as ftme
from fess.builder.samplingStatisticsNew2 import OutfileParser
import time
from multiprocessing.dummy import Pool, TimeoutError #We like a pool of threads
from itertools import islice
import random
import fess.builder.energy as fbe
def get_parser():
    """
    Here all commandline & help-messages arguments are defined.

    :returns: an instance of argparse.ArgumentParser
    """
    parser = argparse.ArgumentParser()
    #Argument(s)
    parser.add_argument('ernwin_dir', type = str, help='A directory containing the step*.coord files.')
    #Options
    parser.add_argument('--fair-ensemble', action='store', help='A folder containing the '
                            'build*.coord files representing a fairly built ensemble. '
                            'Or several ":" seperated folders', type=str)
    parser.add_argument('--reference', type=str, help="The true crystal structure")
    parser.add_argument('--subsample', type=int, default = 1, help="Use every n structures")
    parser.add_argument('--skip-burnin', type=int, default = 0, help="Skip n structures")
    parser.add_argument('--max-bg-samples', type=int, default = 0, help="Use only up to n background structures")
    parser.add_argument('--full-rmsd-matrix', action="store_true", help="Calculate complete RMSD matrix (takes quite long)")
    parser.add_argument('--bg-energy', action="store_true", help="Calculate the energy for the background.")
    parser.add_argument('--ml', type=str, help="Two ,-seperated ml names")
    parser.add_argument('--write-stats-csv', type=str, help="Write this file")
    return parser

def read_cg_traj(fn):
    step = int(re.search(r'step(\d*).coord', fn).group(1))
    cg = ftmc.CoarseGrainRNA(fn)
    return step, cg
def read_cg_bg(fn):
    return fn, ftmc.CoarseGrainRNA(fn)


def sample_background_rmsd(reference):
    oldav = float('inf')
    samples = [random.choice(reference).coords.rmsd_to(random.choice(reference).coords)]
    newav = np.mean(samples)
    while abs(oldav-newav)>0.5:
        oldav = newav
        for i in range(20):
            samples.append(random.choice(reference).coords.rmsd_to(random.choice(reference).coords))
        newav = np.mean(samples)
    return newav

DEFAULT_ENERGY_PREFACTOR=30


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



parser = get_parser()
if __name__=="__main__":
    with open("test", "w") as f:
        f.write(str(1234.2))

    args = parser.parse_args()
    
    #Read the out.log file
    #data = OutfileParser.parse(op.join(args.ernwin_dir,"out.log"))

    # Read all structures
    if args.reference:
        ref = ftmc.CoarseGrainRNA(args.reference)
    else:
        ref = None
    print("Reading traj") 
    traj = {}
    pool = Pool(processes=2)
    file_iter = islice(sorted(glob.iglob(op.join(args.ernwin_dir,"step*.coord"))), args.skip_burnin, None, args.subsample)
        
    for step, cg in pool.imap_unordered(read_cg_traj, file_iter):
        traj[step]=cg
            
    trajectory = ftme.Ensemble(traj, ref)
    print("traj read")
    if args.bg_energy:
        energy_function = fbe.CombinedEnergy([], getDefaultEnergies(trajectory.at_timestep(0)))
    
    sm = lambda : None #Dummy object http://stackoverflow.com/a/2827734/5069869
    ftraj = []
    f_energies = []
    f_fns = []
    if args.fair_ensemble:
        fair_es = args.fair_ensemble.split(":")
        print("Reading BG")               
        for fe in fair_es:
            file_iter = glob.iglob(op.join(fe,"build*.coord"))
            if args.max_bg_samples:
                file_iter = islice(file_iter, 0, args.max_bg_samples)
            for fn, cg in pool.imap(read_cg_bg, file_iter):
                ftraj.append(cg)
                if args.bg_energy:
                    sm.bg = cg
                    e = energy_function.eval_energy(sm, background = False)
                    f_energies.append(e)
                    f_fns.append(fn)
        print("Fair ensemble with {} builds".format(len(ftraj)))        
        if args.bg_energy:
            with open("energies_of_bg", "w") as f:
                for i,e in enumerate(f_energies):
                    f.write(str("{} {}\n".format(f_fns[i], e)))

    else:
        ftraj = None
    
    
    pool.close()
    
    
    if args.write_stats_csv:
        trajectory.create_element_csv(args.write_stats_csv)
        if ftraj:
            bg_e = ftme.Ensemble(ftraj)
            bg_e.create_element_csv(args.write_stats_csv+".bg.csv")
        #sys.exit()
 
    if args.ml:
        for ml in args.ml.split(":"):
            print(args.ml)
            elem1, elem2 = ml.split(",")

            bins = trajectory.view_2d_hist(ftraj, "stat_angle_{}".format(elem1), "stat_angle_{}".format(elem2))
            trajectory.color_by_energy(bins, ftraj, f_energies, "stat_angle_{}".format(elem1), "stat_angle_{}".format(elem2))
            trajectory.view_2d_projection(ftraj, "stat_angle_{}".format(elem1), "stat_angle_{}".format(elem2), cluster=args.full_rmsd_matrix)    
        
            bins = trajectory.view_2d_hist(ftraj, "stat_angle_{}".format(elem1), "cg_distance_{}".format(elem1))
            trajectory.color_by_energy(bins, ftraj, f_energies, "stat_angle_{}".format(elem1), "cg_distance_{}".format(elem1))
            trajectory.view_2d_projection(ftraj, "stat_angle_{}".format(elem1), "cg_distance_{}".format(elem1), cluster=args.full_rmsd_matrix)    
           

            bins = trajectory.view_2d_hist(ftraj, "stat_angle_{}".format(elem2), "cg_distance_{}".format(elem2))
            trajectory.color_by_energy(bins, ftraj, f_energies, "stat_angle_{}".format(elem2), "cg_distance_{}".format(elem2))
            trajectory.view_2d_projection(ftraj, "stat_angle_{}".format(elem2), "cg_distance_{}".format(elem2), cluster=args.full_rmsd_matrix)    
 
            bins = trajectory.view_2d_hist(ftraj, "cg_distance_{}".format(elem1), "cg_distance_{}".format(elem2))
            trajectory.color_by_energy(bins, ftraj, f_energies, "cg_distance_{}".format(elem1), "cg_distance_{}".format(elem2))
            trajectory.view_2d_projection(ftraj, "cg_distance_{}".format(elem1), "cg_distance_{}".format(elem2), cluster=args.full_rmsd_matrix)    
        
            s = "cg_dist_sum_{}_{}".format(elem1, elem2)
            d = "cg_dist_difference_{}_{}".format(elem1, elem2)
            bins = trajectory.view_2d_hist(ftraj, s, d)
            trajectory.color_by_energy(bins, ftraj, f_energies, s, d )
            trajectory.view_2d_projection(ftraj, s, d, cluster=args.full_rmsd_matrix)    
    print("PCA")
    trajectory.ensemble_pca(ftraj)
    trajectory.ensemble_pca(ftraj, False)

    
    print(time.time(),"rmsd - rmsd")    
    bins = trajectory.view_2d_hist(ftraj)
    trajectory.color_by_energy(bins, ftraj, f_energies)
    trajectory.view_2d_projection(ftraj, cluster=args.full_rmsd_matrix)    

    print(time.time(),"rog - rmsd")
    bins = trajectory.view_2d_hist(ftraj, "rog", "rmsd_to_reference")
    trajectory.color_by_energy(bins, ftraj, f_energies, "rog", "rmsd_to_reference")
    trajectory.view_2d_projection(ftraj, "rog", "rmsd_to_reference", cluster=args.full_rmsd_matrix)
    trajectory.view_2d_projection(ftraj, "rog", "rmsd_to_reference", cluster=args.full_rmsd_matrix, circular = True)
    
    print(time.time(),"rog - anisotropy")
    bins = trajectory.view_2d_hist(ftraj, "rog", "anisotropy")
    trajectory.color_by_energy(bins, ftraj, f_energies, "rog", "anisotropy")
    trajectory.view_2d_projection(ftraj, "rog", "anisotropy", cluster=args.full_rmsd_matrix)
    
    if args.full_rmsd_matrix:    
        print(time.time(),"Delta RMSD")
        trajectory.view_delta_rmsd_vs_steps()
        print(time.time(),"Embedding")
        trajectory.view_2d_embedding()
    if ftraj:
        print("The average RMSD between two structures in the background is roughly {:.0f}".format(sample_background_rmsd(ftraj)))
    if ref and ftraj:
        print("The average RMSD between the reference structure and the background is {:.0f}".format(np.mean(ftme.calculate_descriptor_for("rmsd_to_reference", ftraj, ref))))


    if ftraj:
        num_stru=[]
        av_rmsd=[]
        av_true_rmsd=[]
        min_true_rmsd = []
        for i in range(1000, len(ftraj), 500):
            num_stru.append(i)
            av_rmsd.append(sample_background_rmsd(ftraj[:i]))
            rmsds = ftme.calculate_descriptor_for("rmsd_to_reference", ftraj[:i], ref)
            av_true_rmsd.append(np.mean(rmsds))
            min_true_rmsd.append(np.min(rmsds))

        import matplotlib.pyplot as plt
        
        plt.plot(num_stru, av_rmsd, "o-", label="estim. pairwise RMSD")
        plt.plot(num_stru, av_true_rmsd, "o-", label="av. RMSD to target")
        plt.plot(num_stru, min_true_rmsd, "o-", label="min RMSD to target")

        plt.legend(loc = "center right")
        figname = "RMSD_vs_samplesize_{}_background.svg".format(trajectory.at_timestep(0).name)
        plt.savefig(figname)
        log.info("Figure {} created".format(figname))
        plt.clf()
        plt.close()

    num_stru=[]
    av_rmsd=[]
    av_true_rmsd=[]
    min_true_rmsd = []
    for i in range(1000, len(traj), 500):
        num_stru.append(i)
        av_rmsd.append(sample_background_rmsd(trajectory.at_timestep((0,i))))
        rmsds = ftme.calculate_descriptor_for("rmsd_to_reference", trajectory.at_timestep((0,i)), ref)
        av_true_rmsd.append(np.mean(rmsds))
        min_true_rmsd.append(np.min(rmsds))

    import matplotlib.pyplot as plt
    plt.plot(num_stru, av_rmsd, "o-", label="estim. pairwise RMSD")
    plt.plot(num_stru, av_true_rmsd, "o-", label="av. RMSD to target")
    plt.plot(num_stru, min_true_rmsd, "o-", label="min RMSD to target")

    plt.legend(loc = "center right")
    figname = "RMSD_vs_samplesize_{}.svg".format(trajectory.at_timestep(0).name)
    plt.savefig(figname)
    log.info("Figure {} created".format(figname))
    plt.clf()
    plt.close()

        
        
    #Energies only for the steps stored
   # energies = np.array([data["Sampling_Energy"][i] for i in sorted(nr_to_step.values())])
    

