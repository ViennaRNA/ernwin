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
import itertools as it
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.ensemble as ftme
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
    traj = {}
    for fn in glob.iglob(op.join(args.ernwin_dir,"step*.coord")):
        step = int(re.search(r'step(\d*).coord', fn).group(1))
        traj[step] = ftmc.CoarseGrainRNA(fn)
    if args.reference:
        ref = ftmc.CoarseGrainRNA(args.reference)
    else:
        ref = None  
    trajectory = ftme.Ensemble(traj, ref)
    trajectory.view_db_clustering()
    
    
    #Energies only for the steps stored
   # energies = np.array([data["Sampling_Energy"][i] for i in sorted(nr_to_step.values())])
    

