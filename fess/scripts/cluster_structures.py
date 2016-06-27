#!/usr/bin/python
from __future__ import print_function, division

import collections as col
import itertools as it
import numpy as np
import matplotlib.pyplot as plt
import random
import scipy.spatial.distance as ssd
import scipy.cluster.hierarchy as sch
import sys
from optparse import OptionParser
import forgi.threedee.utilities.rmsd as ftur
import forgi.threedee.model.coarse_grain as ftmc
import forgi.utilities.debug as fud
import forgi.threedee.utilities.graph_pdb as cgg
import forgi.threedee.utilities.rmsd as ftur
import forgi.threedee.utilities.vector as ftuv

import sklearn.cluster

#import Pycluster as pc

import scipy.cluster.vq as scv

def cluster_hierarchical(coords, matrix=None, use_heuristic=True):
    # scipy hierarchical is explained in: 
    # https://joernhees.de/blog/2015/08/26/scipy-hierarchical-clustering-and-dendrogram-tutorial/
    dists = np.zeros((len(coords), len(coords)))
    rmsds = np.zeros((len(coords), len(coords)))
    rogs= np.zeros((len(coords)))
    maxgoodrog=float("inf")
    calc=0
    heur=0
    HEUR_FILENUM=50 #If more then this number of files are present and use_heuristic is true, then use a heuristic
    MAX_RMSD=10 #Max RMSD to be considered in a cluster.
    MAX_DRMSD=10 #Max DRMSD to be considered in a cluster.
    print("{} Files".format(len(coords)))
    if use_heuristic and len(coords)>HEUR_FILENUM: #Need to use ROG as heuristic
        for i, c in enumerate(coords):
            rogs[i]=ftur.radius_of_gyration(c)
        for i,j in it.combinations(range(HEUR_FILENUM), r=2):
            dists[i,j] = ftur.drmsd(coords[i], coords[j])
            dists[j,i] = dists[i,j]
            rmsds[i,j] = ftur.rmsd(coords[i], coords[j])
            rmsds[j,i] = rmsds[i,j]
        goodrogs=[abs(rogs[i]-rogs[j]) for  i,j in it.combinations(range(HEUR_FILENUM), 2) if not (dists[i,j]>MAX_DRMSD or rmsds[i,j]>MAX_RMSD)]
        if len(goodrogs)<len(rogs): #We can estimate our heuristic cutoff!
            maxgoodrog=max(abs(rogs[i]-rogs[j]) for  i,j in it.combinations(range(HEUR_FILENUM), 2) if not (dists[i,j]>MAX_DRMSD or rmsds[i,j]>MAX_RMSD))
    for i,j in it.combinations(range(len(coords)), r=2):
        if dists[i,j]>0:
            continue #Already calculated when estimating the heuristic parameter
        if abs(rogs[i]-rogs[j])>maxgoodrog:
            dists[i,j] = 2*MAX_DRMSD
            dists[j,i] = 2*MAX_DRMSD
            rmsds[i,j] = 2*MAX_RMSD
            rmsds[j,i] = 2*MAX_RMSD
            heur+=1
        else:
            dists[i,j] = ftur.drmsd(coords[i], coords[j])
            dists[j,i] = dists[i,j]
            rmsds[i,j] = ftur.rmsd(coords[i], coords[j])
            rmsds[j,i] = rmsds[i,j]
            calc+=1
    fig, ax = plt.subplots(2,2)
    flatdists=np.array([dists[i,j] for i,j in it.combinations(range(len(coords)), 2)])
    flatrogs=np.array([ abs(rogs[i]-rogs[j]) for  i,j in it.combinations(range(len(coords)), 2)])
    flatrmsds=np.array([ abs(rmsds[i,j]) for  i,j in it.combinations(range(len(coords)), 2)])
    cutdists=np.array([dists[i,j] for i,j in it.combinations(range(len(coords)), 2) if dists[i,j]>MAX_DRMSD and rmsds[i,j]>MAX_RMSD])
    cutrogs=np.array([ abs(rogs[i]-rogs[j]) for  i,j in it.combinations(range(len(coords)), 2) if dists[i,j]>MAX_DRMSD and rmsds[i,j]>MAX_RMSD])
    cutrmsds=np.array([ abs(rmsds[i,j]) for  i,j in it.combinations(range(len(coords)), 2) if dists[i,j]>MAX_DRMSD and rmsds[i,j]>MAX_RMSD])
    heurdists=np.array([dists[i,j] for i,j in it.combinations(range(len(coords)), 2) if abs(rogs[i]-rogs[j])>maxgoodrog ])
    heurrogs=np.array([ abs(rogs[i]-rogs[j]) for  i,j in it.combinations(range(len(coords)), 2) if abs(rogs[i]-rogs[j])>maxgoodrog ])
    heurrmsds=np.array([ abs(rmsds[i,j]) for  i,j in it.combinations(range(len(coords)), 2) if abs(rogs[i]-rogs[j])>maxgoodrog ])

    ax[0,0].plot(flatdists, flatrogs, "o")
    ax[1,0].plot(flatrmsds, flatrogs, "o")
    ax[0,1].plot(flatdists, flatrmsds, "o")
    ax[0,0].plot(cutdists, cutrogs, "o")
    ax[1,0].plot(cutrmsds, cutrogs, "o")
    ax[0,1].plot(cutdists, cutrmsds, "o")
    ax[0,0].plot(heurdists, heurrogs, "o")
    ax[1,0].plot(heurrmsds, heurrogs, "o")
    ax[0,1].plot(heurdists, heurrmsds, "o")
    ax[0,0].set_xlabel("drmsd")
    ax[0,0].set_ylabel("delta rog")
    ax[0,1].set_xlabel("drmsd")
    ax[0,1].set_ylabel("RMSDS")
    ax[1,0].set_xlabel("RMSDS")
    ax[1,0].set_ylabel("delta rog")
    plt.show()
    print ("Heuristic: {}, Calculated: {}, Percent Heur: {}".format(heur, calc, heur/(heur+calc)))
    z = sch.linkage(dists, method='complete')
    #print(dists)
    #print (z)
    fig, ax = plt.subplots()
    ax.set_title('Hierarchical Clustering Dendrogram')
    ax.set_xlabel('sample index')
    ax.set_ylabel('distance')
    sch.dendrogram(
        z,
        leaf_rotation=90.,  # rotates the x axis labels
        leaf_font_size=8.,  # font size for the x axis labels
        ax = ax,
        color_threshold = 10,
    )
    plt.axhline(y=10, c='k')

    #fud.pv('dists')
    if matrix is not None:
        np.savetxt(matrix, dists, delimiter=' ', fmt="%.3f")

    ax.set_yscale("symlog", nonposx='clip')
    plt.show()
    flat_clust = sch.fcluster(z, MAX_DRMSD, "distance")
    print(len(set(flat_clust)), "flat clusters")
    print(flat_clust)
    return flat_clust, dists
    
    
    

def distances(s):
    '''
    Compute the distance array for a shape s.
    '''
    ds = [ftuv.vec_distance(p1, p2) for p1,p2 in it.combinations(s, r=2)]
    
    return np.array(ds)

def cluster_kmeans_core(coords, num_clusters=8):
    new_coords = []
    for c in coords:
        dists = distances(c)
        new_coords += [dists]

    return pc.kcluster(new_coords, num_clusters)

def cluster_kmeans(coords, names, topn=0, num_clusters=8):
    labels, error, nfound = cluster_kmeans_core(coords, num_clusters)
    c = col.Counter(labels)

    sorted_labels = sorted(c.keys(), key=lambda x: -c[x])
    #print c

    for l, n in zip(labels, names):
        if l == sorted_labels[topn]:
            print (n)

def cg_fns_to_coords(args):
    structs = []
    coords = []

    for arg in args:
        structs += [ftmc.CoarseGrainRNA(arg)]
        coords += [cgg.bg_virtual_residues(structs[-1])]
    return coords, structs


def main():
    usage = """
    python cluster_structures file1.cg file2.cg ....
    """
    num_args= 1
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')
    parser.add_option('', '--hierarchical', dest='hierarchical', default=False, action='store_true', help='Use hierarchical clustering')
    parser.add_option('', '--topn', dest='topn', default=0, help="Print the entries in the top n'th cluster", type='int')
    parser.add_option('', '--matrix', dest='matrix', default=None, help='Print out the distance matrix', type='str')
    parser.add_option('', '--args', dest='args', default=None, help='Arguments filename', type='str')
    parser.add_option('-n', '--num-structs', dest='num_structs', default=None, help='The number of structures to use', type=int)
    parser.add_option('', '--num-clusters', dest='num_clusters', default=8, help='The number of clusters to use for k-means', type=int)

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    random.shuffle(args)

    if options.num_structs is not None:
        args = args[:options.num_structs]

    coords, structs = cg_fns_to_coords(args)

    if options.args is not None:
        with open(options.args, 'w') as f:
            f.write(" ".join(args))
    
    if options.hierarchical:
        clusters, dists = cluster_hierarchical(coords, options.matrix)
    else:
        cluster_kmeans(coords, args, options.topn, 
                       num_clusters=options.num_clusters)
    clustered_filenames=col.defaultdict(list)
    for i, filename in enumerate(args):
        clustered_filenames[clusters[i]].append((filename,i))
    for key in sorted(clustered_filenames.keys(), key=lambda x: -len(clustered_filenames[x])):
        print ("CLUSTER {} with {} members".format(key, len(clustered_filenames[key])))
        bestscore=float("inf")
        bestfn = None
        for i, fn in enumerate(clustered_filenames[key]):
            print("\t", fn[0])
            score=0
            for j, fn2 in enumerate(clustered_filenames[key]):
                if i!=j:
                    score+=dists[fn[1],fn2[1]]
            if score<bestscore:
                bestscore=score
                bestfn=fn[0]
        print("\t REPRESENTATIVE: ", bestfn)
        
        
            

if __name__ == '__main__':
    main()

