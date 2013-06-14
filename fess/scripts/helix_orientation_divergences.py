#!/usr/bin/python

import collections as c
import scipy.cluster.hierarchy as sch
import itertools as it
import random
import numpy as np
import math as m
import corgy.builder.stats as cbs

from matplotlib import cm
from matplotlib.ticker import LinearLocator
import sys, pandas as pa
import scipy.stats as ss
from optparse import OptionParser
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d, Axes3D

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

import corgy.utilities.debug as cud
import corgy.utilities.vector as cuv

def get_certain_angle_stats(stats, angle_type):
    '''
    Return an index containing only the statisitcs corresponding to
    this type of angle.
    '''
    ss = angle_type
    ss_r = stats[stats["s1"] == ss[0]]
    ss_r = ss_r[ss_r["s2"] == ss[1]]
    ss_r = ss_r[ss_r["atype"] == ss[2]]
    return ss_r

def spherical_distance(phi1, phi2, theta1, theta2):
    '''
    Returns the distance between the coordinates on the unit sphere.
    '''
    cos = (math.sin(phi1)*math.sin(phi2)*math.cos(theta1 - theta2) + 
           math.cos(phi1)*math.cos(phi2))
    arc = math.acos( cos )

    # Remember to multiply arc by the radius of the earth 
    # in your favorite set of units to get length.
    return arc

def get_nearest_dimension_sizes(bulge_size, stat_counts, min_entries = 20):
    '''
    Classify all the bulge sizes according to their distance from the
    one being analyzed. Then pick the next n closest such that the total
    number of statistics for all of the selected bulge sizes exceeds
    min_entries.
    '''
    sizes = []
    for key in stat_counts.keys():
        if key[2] != bulge_size[2]:
            # not the same angle type
            continue
        else:
            dist = (key[1] - bulge_size[1]) ** 2 + (key[0] - bulge_size[0]) ** 2 
            sizes += [tuple([dist] + list(key))]

    sizes.sort()
    selected_sizes = []
    cumulative_num_entries = 0

    for s in sizes:
        cumulative_num_entries += stat_counts[tuple(s[1:])]
        selected_sizes += [tuple(s[1:])]

        if cumulative_num_entries >= min_entries:
            break

    return (selected_sizes, cumulative_num_entries)

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
    ./helix_orienation_divergences.py

    Analyze how much the helix-helix orientations diverge between two data sets.
    """
    num_args=0
    parser = OptionParser()

    parser.add_option('-r', '--resolution', dest='resolution', default=10, help="The resolution of the resulting plot", type='int')
    parser.add_option('-a', '--angle', dest='angle', default=0, help="The angle of the camera", type='float')
    parser.add_option('-f', '--fig-name', dest='fig_name', default='', help="The name of the file to save the figure to. If it is not specified, the figure will not be saved", type='str')
    parser.add_option('-i', '--interior_loops', dest='interior_loops', default=False, help='Cluster only the interior loops', action='store_true')
    parser.add_option('-m', '--multi_loops', dest='multi_loops', default=False, help='Cluster only the interior loops', action='store_true')

    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    column_names = ['type', 'pdb', 's1', 's2', 'u', 'v', 't', 'r', 'u1', 'v1', 'atype', 'something1', 'something2', 'sth3', 'sth4']
    real_stats = pa.read_csv('fess/stats/temp.real.stats', header=None, sep=' ', names=column_names, engine='python')
    sampled_stats = pa.read_csv('fess/stats/temp.sampled.stats', header=None, sep=' ', names=column_names, engine='python')

    real_stats = real_stats[real_stats["type"] == "angle"]
    real_stat_dims = map(tuple, real_stats[['s1', 's2', 'atype']].as_matrix())

    # count how many statistics we have for each statistic type
    stat_counts = c.defaultdict(int)
    for sc in real_stat_dims:
        stat_counts[sc] += 1

    cud.pv('stat_counts')
    histograms = dict()
    for b in stat_counts.keys():
        if b[2] != 2.:
            # only look at type 2 angles
            continue

        if options.interior_loops:
            if b[0] == 1000 or b[1] == 1000:
                continue
        if options.multi_loops:
            if b[0] != 1000 and b[1] != 1000:
                continue

        (selected_sizes, count) = get_nearest_dimension_sizes(b, stat_counts, 1)

        if count < 3:
            continue

        cud.pv('b, selected_sizes')

        combined_real = []

        # get the statistics that correspond to the selected sampled sizes
        for ss in selected_sizes:
            ss_r = get_certain_angle_stats(real_stats, ss)

            combined_real += list(ss_r[['u','v']].as_matrix())

        num_points = len(combined_real)
        combined_real = np.array(combined_real)
        #histograms[b] = (np.histogram2d(combined_real[:,0], combined_real[:,1], range=[[0, m.pi], [-m.pi, m.pi]])[0] + 0.5) / float(num_points)
        histograms[b] = combined_real

    dists = []
    named_dists = dict()
    pp_dists = dict()
    for k1, k2 in it.combinations(histograms.keys(), 2):
        per_point_distances = []
        for p1 in histograms[k1]:
            point_distances = []
            for p2 in histograms[k2]:
                point_distances += [cuv.magnitude(p1 - p2)]
            per_point_distances += [min(point_distances)]

        for p2 in histograms[k2]:
            point_distances = []
            for p1 in histograms[k1]:
                point_distances += [cuv.magnitude(p1-p2)]
            per_point_distances += [min(point_distances)]

        dists += [max(per_point_distances)]
        named_dists[(k1,k2)] = max(per_point_distances)
        pp_dists[(k1,k2)] = per_point_distances

        '''
        kl = histograms[k1] * (histograms[k1] / histograms[k2])
        kl = sum(map(sum, kl))
        dists += [kl]
        '''

    cud.pv('dists')
    Z = sch.complete(dists)
    cud.pv('Z')
    sch.dendrogram(Z, labels = histograms.keys(), leaf_rotation=90)
    plt.subplots_adjust(bottom=0.25)
    
    plt.show()
    sys.exit(1)

    k1 = (6,7,2)
    k2 = (5,6,2)

    rs = get_certain_angle_stats(real_stats, k1)
    ss = get_certain_angle_stats(real_stats, k2)

    cud.pv('named_dists[(k1,k2)]')
    cud.pv('pp_dists[(k1,k2)]')

    real_us = rs[['u', 'v']].as_matrix()
    sampled_us = ss[['u','v']].as_matrix()

    U_r = real_us[:,0]
    V_r = real_us[:,1]

    U_s = sampled_us[:,0]
    V_s = sampled_us[:,1]

    total_r = len(U_r)
    total_s = len(U_s)

    hr = np.histogram2d(U_r, V_r)
    hs = np.histogram2d(U_s, V_s)

    pseudo_r = (hr[0] + 1) / total_r
    pseudo_s = (hs[0] + 1) / total_r
    kl = pseudo_r * (pseudo_r / pseudo_s)
    cud.pv('kl')
    cud.pv('sum(map(sum, kl))')

    X_r = np.sin(U_r) * np.cos(V_r)
    Y_r = np.sin(U_r) * np.sin(V_r)
    Z_r = np.cos(U_r)

    r = 1.
    X_s = r * np.sin(U_s) * np.cos(V_s)
    Y_s = r * np.sin(U_s) * np.sin(V_s)
    Z_s = r * np.cos(U_s)

    cud.pv('real_us')

    real_us_orig = np.copy(real_us)
    sampled_us_orig = np.copy(sampled_us)

    print len(real_us), len(sampled_us)

    fig = plt.figure(figsize=(10,10))
    ax = Axes3D(fig)

    a = Arrow3D([-1.3,1.3],[0,0],[0,0], mutation_scale=20, lw=5, arrowstyle="-|>", color="g")
    ax.add_artist(a)

    ax.plot(X_r, Y_r, Z_r, 'bo', alpha=0.3)
    ax.plot(X_s, Y_s, Z_s, 'ro', alpha=0.3)

    u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
    x=np.cos(u)*np.sin(v)
    y=np.sin(u)*np.sin(v)
    z=np.cos(v)
    ax.plot_wireframe(x, y, z, color="y")

    #surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=colors,
    #       linewidth=0, antialiased=False)

    ax._axis3don=False
    ax.set_zlim3d(-1, 1)
    ax.w_zaxis.set_major_locator(LinearLocator(6))
    ax.view_init(0, options.angle)

    '''
    plt.subplots_adjust(left=0.4, right=0.9, top=0.9, bottom=0.1)

    for i in xrange(0, 360, 40):
        savefig("fig%d.png", (i))
    '''

    '''
    sm = cm.ScalarMappable(cmap=cm.jet)
    sm.set_array(W)
    fig.colorbar(sm)
    '''

    if options.fig_name != "":
        plt.savefig(options.fig_name, bbox_inches='tight')
    else:
        plt.show()

if __name__ == '__main__':
    main()

