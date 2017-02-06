#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip) #future package
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)
import logging    
logging.basicConfig(level=logging.INFO)
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
                            'build*.coord files representing a fairly built ensemble.', type=str)
    parser.add_argument('--reference', type=str, help="The true crystal structure")
    parser.add_argument('--subsample', type=int, default = 1, help="Use every n structures")

    return parser

def read_cg_traj(fn):
    step = int(re.search(r'step(\d*).coord', fn).group(1))
    cg = ftmc.CoarseGrainRNA(fn)
    return step, cg
def read_cg_bg(fn):
    return ftmc.CoarseGrainRNA(fn)


def sample_background_rmsd(reference):
    oldav = float('inf')
    samples = [random.choice(reference).coords.rmsd_to(random.choice(reference).coords)]
    newav = np.mean(samples)
    while abs(oldav-newav)>2:
        oldav = newav
        for i in range(20):
            samples.append(random.choice(reference).coords.rmsd_to(random.choice(reference).coords))
        newav = np.mean(samples)
    return newav
    
parser = get_parser()
if __name__=="__main__":

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
    file_iter = islice(sorted(glob.iglob(op.join(args.ernwin_dir,"step*.coord"))), 0, None, args.subsample)
        
    for step, cg in pool.imap_unordered(read_cg_traj, file_iter):
        traj[step]=cg
            
    trajectory = ftme.Ensemble(traj, ref)
    del traj
    print("traj read")
    if args.fair_ensemble:
        print("Reading BG")        
        ftraj = []
        file_iter = glob.iglob(op.join(args.fair_ensemble,"build*.coord"))
        for cg in pool.imap_unordered(read_cg_bg, file_iter):
            ftraj.append(cg) #Arbitrary order
        print("Fair ensemble with {} builds".format(len(ftraj)))
    else:
        ftraj = None
        
    pool.close()
    
    print(time.time(),"rmsd - rmsd")
    trajectory.view_db_clustering(ftraj)
    print(time.time(),"rog - rmsd")
    trajectory.view_db_clustering(ftraj, "rog", "rmsd_to_reference")
    print(time.time(),"rog - anisotropy")
    trajectory.view_db_clustering(ftraj, "rog", "anisotropy")
    print(time.time(),"Delta RMSD")
    trajectory.view_delta_rmsd_vs_steps()
    print(time.time(),"Embedding")
    trajectory.view_2d_embedding()
    if ftraj:
        print("The average RMSD between two structures in the background is roughly {:.0f}".format(sample_background_rmsd(ftraj)))
    if ref and ftraj:
        print("The average RMSD between the reference structure and the background is {:.0f}".format(np.mean(calculate_descriptor_for("rmsd_to_reference", ftraj, ref))))

    #Energies only for the steps stored
   # energies = np.array([data["Sampling_Energy"][i] for i in sorted(nr_to_step.values())])
    

