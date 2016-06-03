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
import forgi.threedee.utilities.rmsd as ftur
import forgi.threedee.utilities.graph_pdb as ftug
from fess import data_file
import numpy as np
import scipy.cluster.hierarchy as hcluster
import matplotlib.pyplot as plt
import copy
import collections as col
import sklearn.cluster

def generateParser():
    parser=argparse.ArgumentParser( description="Estimate the size of the solution space "
                                                "for the sampling of a coarse grained structure.")
    parser.add_argument("cg", type=str, action="store", help="A cg file")
    parser.add_argument("--stats-file", type=str, help="What stats file to use", default=data_file("stats/all.stats"))
    parser.add_argument("-t", "--tolerance", type=float, help="Element-wise relative tolerance of r,u, v and t for considering two stats as equal.", default=0.1)
    parser.add_argument("-r", "--rmsd", action="store_true", help="Use RMSD of whole structure to identify similar stats for individual coarse grain elements (takes more time)")
    parser.add_argument("-c", "--cluster", action="store_true", help="Calculate pairwise distances of stats for clustering")
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

class Comparison(object):
    def __init__(self, reference_cg, tolerance=0.1, rmsd=False):
        self.tol=tolerance
        self.rmsd=rmsd
        self.restore_cg(reference_cg)
    def restore_cg(self, reference_cg):
        self.ref_sm = fbm.SpatialModel(copy.deepcopy(reference_cg))
        self.ref_sm.load_sampled_elems()
        self.ref_vress = self.ref_sm.bg.get_ordered_virtual_residue_poss()
    def compare(self, d, other_stats):
        if not self.rmsd:
            stats = self.ref_sm.elem_defs[d]
            return is_similar(stats, other_stats, self.tol)
        else:
            self.ref_sm.elem_defs[d]=other_stats
            self.ref_sm.traverse_and_build(start=d)
            curr_vress=self.ref_sm.bg.get_ordered_virtual_residue_poss()
            rmsd=ftur.rmsd(self.ref_vress, curr_vress)
            return rmsd<self.tol

parser=generateParser()
if __name__=="__main__":
    args = parser.parse_args()
    cg=ftmc.CoarseGrainRNA(args.cg)
    sm=fbm.SpatialModel(cg, ftms.get_conformation_stats(args.stats_file))
    sm.load_sampled_elems()
    maxNum=0
    maxLoop=0
    sumStats=0
    comp = Comparison(sm.bg, args.tolerance, args.rmsd)
    for define in cg.defines:
        similar_stats=0
        allstats = sm.conf_stats.sample_stats(cg, define)
        numStats = len(allstats)
        if define[0] in "mi" and numStats:        
            if args.cluster:
                distance=np.full((len(allstats),len(allstats)),np.NAN)
            for i, stat in enumerate(allstats):
                if args.cluster:
                    for j, stat2 in enumerate(allstats):
                       distance[i,j]=stat.diff(stat2) 
                if comp.compare(define, stat):
                    similar_stats+=1
            if args.cluster:
                np.savetxt(define+'.distances.out', distance, delimiter='\t', header="\t".join(s.pdb_name+",{},{}".format(s.r1, s.u1) for s in allstats ))
                db = sklearn.cluster.DBSCAN(metric="precomputed", eps=2.1, min_samples=2).fit(distance)
                labels = db.labels_
                n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
                n_outliers = len([x for x in labels if x==-1])
                n_largest_clust = col.Counter(labels).most_common(1)
                print ("DBSCAN clustering: {} clusters, {} outliers, {} in largest cluster".format(n_clusters, n_outliers, n_largest_clust[0][1]))
            effectiveNum = numStats/similar_stats
            print ("{}: {} different stats, 1/{} success rate".format(define, numStats, effectiveNum))
            comp.restore_cg(sm.bg)
        else:
            effectiveNum = numStats
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

