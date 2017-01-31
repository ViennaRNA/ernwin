#!/usr/bin/python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip) #future package
from future.builtins.disabled import (apply, cmp, coerce, execfile,
                             file, long, raw_input, reduce, reload,
                             unicode, xrange, StandardError)

import argparse
import glob
import os.path as op
import re
import numpy as np
import itertools as it
import forgi.threedee.model.coarse_grain as ftmc
from sklearn.cluster import DBSCAN
from fess.builder.samplingStatisticsNew2 import OutfileParser
import matplotlib.pyplot as plt

def distance_copy_previouse_line(distance_matrix, i):
    print ("L", end="")
    distance_matrix[i, :] = distance_matrix[i-1,:]
    distance_matrix[:, i] = distance_matrix[:, i-1]

def distance_copy_previouse_field(distance_matrix, i, j):
    print ("F", end="")
    distance_matrix[i,j] = distance_matrix[j,i] = distance_matrix[i,j-1]

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

    return parser

parser = get_parser()
if __name__=="__main__":
    args = parser.parse_args()
    
    #Read the out.log file
    data = OutfileParser.parse(op.join(args.ernwin_dir,"out.log"))

    # Read all structures
    trajectory = {}
    for fn in glob.iglob(op.join(args.ernwin_dir,"step*.coord")):
        step = int(re.search(r'step(\d*).coord', fn).group(1))
        trajectory[step] = ftmc.CoarseGrainRNA(fn)
    
    nr_to_step = { i: step for i, step in enumerate(sorted(trajectory.keys()))}

    if args.reference:
        ref_cg = ftmc.CoarseGrainRNA(args.reference)
    else:
        ref_cg = trajectory[nr_to_step[0]]
        
    # Create the pairwise distance matrix using cg-RMSD
    distance_matrix = np.zeros((len(trajectory), len(trajectory)))*np.nan
    print("Calculating distance matrix for {} structures".format(len(trajectory)))
    #Calculate RMSD to precursor
    for i in range(1, len(trajectory)):
        j=i-1
        distance_matrix[i,j] = distance_matrix[j,i] = trajectory[nr_to_step[i]].coords.rmsd_to(trajectory[nr_to_step[j]].coords)
    #Calculate rest of distance matrix. If two subsequent entries are the same, copy the values only.
    for i in range(len(trajectory)):
        if i>0 and distance_matrix[i,i-1] < 10**-10: #(Almost) 0 distance
            distance_copy_previouse_line(distance_matrix, i)
            continue
        for j in range(i+2, len(trajectory)):
            if distance_matrix[j,j-1] < 10**-10: #(Almost) 0 distance
                distance_copy_previouse_field(distance_matrix, i, j)
            else:
                distance_matrix[i,j] = distance_matrix[j,i] = trajectory[nr_to_step[i]].coords.rmsd_to(trajectory[nr_to_step[j]].coords)

    print("Saving distance matrix")
    np.savetxt("distance_matrix.txt", distance_matrix)        
    o_rmsd = []
    
    for i in range(len(trajectory)):
        o_rmsd.append(trajectory[nr_to_step[i]].coords.rmsd_to(ref_cg.coords))
    o_rmsd = np.array(o_rmsd)
        
    db = DBSCAN(eps=0.5, min_samples=10, metric="precomputed", n_jobs = 4).fit(distance_matrix)
    
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    
    unique_labels = set(labels)
    colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
    
    #Energies only for the steps stored
    energies = np.array([data["Sampling_Energy"][i] for i in sorted(nr_to_step.values())])
    
    for k, col in zip(unique_labels, colors):
        if k == -1:
            # Black used for noise.
            col = 'k'
        class_member_mask = (labels == k)
        plt.plot( o_rmsd[class_member_mask & core_samples_mask],
                  energies[class_member_mask & core_samples_mask],
                  'o', markerfacecolor=col, markeredgecolor='k', markersize=8
                  )
        plt.plot( o_rmsd[class_member_mask & ~core_samples_mask],
                  energies[class_member_mask & ~core_samples_mask],
                  'o', markerfacecolor=col, markeredgecolor='k', markersize=4
                )
    plt.show()