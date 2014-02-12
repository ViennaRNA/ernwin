#!/usr/bin/python

import itertools as it
import numpy as np
import matplotlib.pyplot as plt
import scipy.spatial.distance as ssd
import scipy.cluster.hierarchy as sch
import sys
from optparse import OptionParser

import forgi.threedee.utilities.rmsd as ftur
import forgi.threedee.model.coarse_grain as ftmc
import forgi.utilities.debug as fud
import forgi.threedee.utilities.graph_pdb as cgg

import scipy.cluster.vq as scv

def main():
    usage = """
    python cluster_structures file1.cg file2.cg ....
    """
    num_args= 1
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    structs = []
    coords = []

    for arg in args:
        structs += [ftmc.CoarseGrainRNA(arg)]
        coords += [cgg.bg_virtual_residues(structs[-1])]

    dists = np.zeros((len(coords), len(coords)))
    for i,j in it.combinations(range(len(coords)), r=2):
        dists[i][j] = ftur.centered_rmsd(coords[i], coords[j])
        dists[j][i] = dists[i][j]

    coords = np.array(coords)
    cl = sch.linkage(dists)

    #fud.pv('dists')
    for i,a in enumerate(args):
        print i, a
    sch.dendrogram(cl)
    fud.pv('cl')
    plt.show()


if __name__ == '__main__':
    main()

