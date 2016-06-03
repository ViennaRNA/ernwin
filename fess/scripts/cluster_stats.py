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
import forgi.threedee.utilities.vector as ftuv
import forgi.threedee.utilities.graph_pdb as ftug
from fess import data_file
import numpy as np
import scipy.cluster.hierarchy as hcluster
import matplotlib.pyplot as plt
import copy, sys
import collections as col
import sklearn.cluster

def generateParser():
    parser=argparse.ArgumentParser( description="Cluster angle stats.")
    parser.add_argument("--stats-file", type=str, help="What stats file to use", default="stats/all_filtered.stats")
    parser.add_argument("-p", "--plot", action="store_true", help="Plot some clusters")
    parser.add_argument("-t", "--threshold", type=float, help="Threshold for the Birch clustering algorithm", default=0.3)
    return parser


parser=generateParser()
if __name__=="__main__":
    args = parser.parse_args()
    #Get the stats object
    conf_stats = ftms.get_conformation_stats(data_file(args.stats_file))
    angle_stats = conf_stats.angle_stats
    colors = plt.cm.Spectral(np.linspace(0, 1, 10))
    colors=np.array(list(colors)*500)
    markers="o"*10+"v"*10+"<"*10+">"*10+"."*10+"s"*10+"^"*10+"8"*10+"*"*10+"+"*10
    markers=markers*500
    for key in angle_stats.keys():
            point_to=[]
            to_clust=[] 
            angles=[]
            for stat in angle_stats[key]:
                stem_start = np.array(ftuv.spherical_polar_to_cartesian([stat.r1, stat.u1, stat.v1]))
                #stem_vec = np.array(ftuv.spherical_polar_to_cartesian([1, stat.u, stat.v]))
                #point1 = stem_start+stem_vec
                stem_vec = np.array(ftuv.spherical_polar_to_cartesian([1, stat.u, stat.v]))
                #point2 = stem_start+stem_vec
                #stem_vec = np.array(ftuv.spherical_polar_to_cartesian([100, stat.u, stat.v]))/100
                #point3 = stem_start+stem_vec
                #point_to.append(list(point1)+list(point2)+list(point3))
                point_to.append(list(stem_start)+list(stem_vec))
                to_clust.append(list(stem_vec))
                angles.append([stat.u, stat.v])

            point_to=np.array(point_to)
            if len(point_to)<10: continue
            #db = sklearn.cluster.DBSCAN(eps=args.threshold, min_samples=3).fit(to_clust)
            #labels_db = db.labels_
            #bandwidth = sklearn.cluster.estimate_bandwidth(point_to, quantile=0.2, n_samples=500)
            #bandwidth=2
            #ms = sklearn.cluster.MeanShift(bandwidth=bandwidth, bin_seeding=True).fit(point_to)
            #cluster_centers = ms.cluster_centers_
            #labels = ms.labels_
            bi = sklearn.cluster.Birch(n_clusters=None, threshold=args.threshold, branching_factor=200, copy=True).fit(to_clust)
            labels = bi.labels_
            n_clusters = len(set(labels)) - (1 if -1 in labels else 0)
            n_outliers = len([x for x in labels if x==-1])
            count_labs = col.Counter(labels)
            n_largest_clust = count_labs.most_common(1)
            n_singletons = len([x for x in count_labs if count_labs[x]==1])
            if len(point_to)>1000:
                   print ("{} - {}: {} clusters. {} in largest cluster. "
                   "{} singletons".format(len(point_to), key, n_clusters, n_largest_clust[0][1],
                   n_singletons), file=sys.stderr)
            
            if args.plot and  n_clusters >0 and n_largest_clust[0][1]>1:
                fig, ax = plt.subplots(2,2)
            for label in sorted(set(labels)):
                if label==-1:
                    print("# OUTLIERS:")
                else:
                    print("# Cluster {} for {}:".format(label, key))
                class_member_mask = (labels == label)
                xy = point_to[class_member_mask] 
                angles=np.array(angles)
                uv = angles[class_member_mask]
                if args.plot and n_largest_clust[0][1]>1 and np.count_nonzero(class_member_mask)>5:
                    ax[0,0].plot(xy[:, 0], xy[:, 1], markers[label], markersize=14, markerfacecolor=colors[label], label=str(label))
                    ax[0,1].plot(uv[:, 0], uv[:, 1], markers[label], markersize=14, markerfacecolor=colors[label], label=str(label))
                    ax[1,0].plot(xy[:, 3], xy[:, 4], markers[label], markersize=14, markerfacecolor=colors[label], label=str(label))
                    ax[1,1].plot(xy[:, 3], xy[:, 5], markers[label], markersize=14, markerfacecolor=colors[label], label=str(label))
                    #ax[2,0].plot(xy[:, 6], xy[:, 7], 'o', markersize=14, markerfacecolor=colors[label], label=str(label))
                    #ax[2,1].plot(xy[:, 6], xy[:, 8], 'o', markersize=14, markerfacecolor=colors[label], label=str(label))
                elif args.plot and n_largest_clust[0][1]>1:
                    ax[0,0].plot(xy[:, 0], xy[:, 1], 'o', markersize=2, color="magenta", mec="magenta", label=str(label))
                    ax[0,1].plot(uv[:, 0], uv[:, 1], 'o', markersize=2, color="magenta", mec="magenta", label=str(label))
                    ax[1,0].plot(xy[:, 3], xy[:, 4], 'o', markersize=2, color="magenta", mec="magenta", label=str(label))
                    ax[1,1].plot(xy[:, 3], xy[:, 5], 'o', markersize=2, color="magenta", mec="magenta", label=str(label))
                    #ax[2,0].plot(xy[:, 6], xy[:, 7], 'o', markersize=2, color="magenta", mec="magenta", label=str(label))
                    #ax[2,1].plot(xy[:, 6], xy[:, 8], 'o', markersize=2, color="magenta", mec="magenta", label=str(label))
                for index in class_member_mask.nonzero()[0]: 
                    print(angle_stats[key][index])
            if args.plot and n_clusters>0 and n_largest_clust[0][1]>1:
                plt.show()
