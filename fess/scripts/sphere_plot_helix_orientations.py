#!/usr/bin/python

import numpy as np
from matplotlib import cm
from matplotlib.ticker import LinearLocator
import sys, pandas as pa
import scipy.stats as ss
from optparse import OptionParser
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    '''
    From:
        
    http://stackoverflow.com/questions/11140163/python-matplotlib-plotting-a-3d-cube-a-sphere-and-a-vector
    '''

    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


def main():
    usage = """
    ./sphere_plot_helix_orientations.py 
    """
    num_args=0
    parser = OptionParser()

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    column_names = ['type', 'pdb', 's1', 's2', 'u', 'v', 't', 'r', 'u1', 'v1', 'atype', 'something1', 'something2', 'sth3', 'sth4']
    real_stats = pa.read_csv('fess/stats/temp.real.stats', header=None, sep=' ', names=column_names, engine='python')
    sampled_stats = pa.read_csv('fess/stats/temp.sampled.stats', header=None, sep=' ', names=column_names, engine='python')

    real_stats = real_stats[real_stats["type"] == "angle"]
    real_us = real_stats[['u', 'v']].as_matrix()
    sampled_us = sampled_stats[['u','v']].as_matrix()

    print "len(real_us)", len(real_us)
    print "len(sampled_us)", len(sampled_us)

    k_r = ss.gaussian_kde(real_us.T)
    k_s = ss.gaussian_kde(sampled_us.T)

    U,V = np.mgrid[0:3.14:10j, -3.14:3.14:10j]
    positions = np.vstack([U.ravel(), V.ravel()])
    vals_r = np.log(k_r(positions))
    vals_s = np.log(k_s(positions))
    vals = vals_r - vals_s

    fig = plt.figure()
    ax = Axes3D(fig)

    ulen = len(U)
    vlen = len(V)

    X = np.sin(U) * np.cos(V)
    Y = np.sin(U) * np.sin(V)
    Z = np.cos(U)


    colors = cm.coolwarm(np.reshape(vals.T, U.shape))

    a = Arrow3D([-1.3,1.3],[0,0],[0,0], mutation_scale=20, lw=2, arrowstyle="-|>", color="g")
    ax.add_artist(a)
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=colors,
            linewidth=0, antialiased=False)


    ax.set_zlim3d(-1, 1)
    ax.w_zaxis.set_major_locator(LinearLocator(6))

    plt.show()

if __name__ == '__main__':
    main()

