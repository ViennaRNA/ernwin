from __future__ import print_function, absolute_import, division
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, open, pow, range, round,
                      str, super, zip) 


import argparse, random, math
import forgi.utilities.stuff as fus
import fess.builder.models as fbm
import fess.builder.energy as fbe
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.stats as ftms
import forgi.threedee.utilities.graph_pdb as ftug
from fess import data_file
import numpy as np
import scipy.cluster.hierarchy as hcluster
import matplotlib.pyplot as plt

def generateParser():
    parser=argparse.ArgumentParser( description="Estimate the size of the solution space "
                                                "for the sampling of a coarse grained structure.")
    parser.add_argument("cg", type=str, action="store", help="A cg file")
    parser.add_argument("--stats-file", type=str, help="What stats file to use", default="stats/all_filtered.stats")
    parser.add_argument("-t", "--tolerance", type=float, help="Element-wise relative tolerance of r,u, v and t for considering two stats as equal.", default=0.1)
    parser.add_argument("-c", "--cluster", action="store_true", help="Cluster the angle stats")
    return parser


def is_similar(stat1, stat2, cutoff=0.1):
    if stat1.dim1 != stat2.dim1 or stat1.dim2 != stat2.dim2:
        return False
    if not np.allclose(stat1.u, stat2.u, rtol=cutoff):
        return False
    if not np.allclose(stat1.v, stat2.v, rtol=cutoff):
        return False
    if not np.allclose(stat1.t, stat2.t, rtol=cutoff):
        return False
    if not np.allclose(stat1.r1, stat2.r1, rtol=cutoff):
        return False
    if not np.allclose(stat1.u1, stat2.u1, rtol=cutoff):
        return False
    if not np.allclose(stat1.v1, stat2.v1, rtol=cutoff):
        return False

    return True
parser=generateParser()
if __name__=="__main__":
    args = parser.parse_args()
    cg=ftmc.CoarseGrainRNA(args.cg)
    sm=fbm.SpatialModel(cg, ftms.get_conformation_stats(data_file(args.stats_file)))
    sm.load_sampled_elems()
    maxNum=0
    maxLoop=0
    sumStats=0

    for define in cg.defines:
        similar_stats=0
        allstats = sm.conf_stats.sample_stats(cg, define)
        numStats = len(allstats)
        if define[0] in "mi" and numStats:
            if args.cluster:
                allstatsvec = np.array([ [ stats.u, stats.v, stats.t, stats.r1, stats.u1, stats.v1 ] 
                                                                         for stats in allstats])
                mean=np.mean(np.abs(allstatsvec), 0)
                allstatsvec=allstatsvec/mean
                thresh = 1
                clusters = hcluster.fclusterdata(allstatsvec, thresh)
                print (np.max(clusters), "Clusters")
            for stat in allstats:
                if is_similar(stat, sm.elem_defs[define], args.tolerance):
                    similar_stats+=1
            effectiveNum = numStats/similar_stats
            print ("{}: {} different stats, 1/{} success rate".format(define, numStats, effectiveNum))
        else:
            effectiveNum == numStats
            print ("{}: {} different stats".format(define, numStats))
        sumStats+=numStats
        if numStats>maxNum:
            maxNum=numStats
        if define[0] in ["m","i"] and effectiveNum>maxLoop:
            maxLoop=effectiveNum
    # The minimal number of sampling steps needed is estimated as follows:
    # The expectation value for the coverage of sampling with replacement is n*(1-(1-1/n)**x) where
    # n is the size of the set drawn from and x is the number of draws.
    # We want at least 50% coverage for the coarse-grained element with the largest number of 
    # possibilities that is a interior or multi loop.
    # We divide this by the probability that this element is changed in a sampling step [1/len(defines)]
    
    needed=math.log(0.5)/math.log(1-1/maxLoop)*len(cg.defines)
    print ("{} cg-elements, {} total stats, {} [effective] for worst angle stats, {} [total] including stems, more than {} sampling steps needed for 50% success chance.".format(
              len(cg.defines), sumStats, maxLoop, maxNum, int(math.ceil(needed))
          ))
    needed90=math.log(0.1)/math.log(1-1/maxLoop)*len(cg.defines)
    print ("{} sampling steps needed for 90% success chance.".format(int(math.ceil(needed90))))

