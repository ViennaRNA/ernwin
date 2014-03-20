#!/usr/bin/python

import collections as col
import itertools as it
import numpy as np
import matplotlib
matplotlib.use("Agg")
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

import Pycluster as pc

import scipy.cluster.vq as scv

def cluster_hierarchical(coords, matrix=None):
    dists = np.zeros((len(coords), len(coords)))
    for i,j in it.combinations(range(len(coords)), r=2):
        dists[i][j] = ftur.centered_drmsd(coords[i], coords[j])
        dists[j][i] = dists[i][j]

    #fud.pv('dists')
    if matrix is not None:
        np.savetxt(matrix, dists, delimiter=' ', fmt="%.3f")

    coords = np.array(coords)
    cl = sch.linkage(dists)

    #fud.pv('dists')
    '''
    for i,a in enumerate(args):
        print i, a
    '''
    sch.dendrogram(cl)
    np.set_printoptions(precision=3, suppress=True)
    #fud.pv('cl')
    plt.show()

    return cl

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
            print n

def cg_fns_to_coords(args):
    structs = []
    coords = []

    for arg in args:
        structs += [ftmc.CoarseGrainRNA(arg)]
        coords += [cgg.bg_virtual_residues(structs[-1])]
    return coords


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

    coords = cg_fns_to_coords(args)

    if options.args is not None:
        with open(options.args, 'w') as f:
            f.write(" ".join(args))
    
    if options.hierarchical:
        cluster_hierarchical(coords, options.matrix)
    else:
        cluster_kmeans(coords, args, options.topn, 
                       num_clusters=options.num_clusters)

if __name__ == '__main__':
    main()

