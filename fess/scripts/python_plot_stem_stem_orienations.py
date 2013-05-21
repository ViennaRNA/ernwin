#!/usr/bin/python

import matplotlib.pyplot as plt
import pandas as pa
import math as m
import sys, numpy as np
import scipy.stats as ss
from optparse import OptionParser

def main():
    usage = """
    ./python_plot_stem_stem_orientations.py
    """
    num_required_args = 0
    parser = OptionParser()

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_required_args:
        parser.print_help()
        sys.exit(1)

    column_names = ['dist', 'out_of_plane', 'in_plane', 'vec_angle', 'lateral_offset', 'ortho_offset']

    real_stats = pa.read_csv('fess/stats/stem_stem_orientations.csv', header=None, sep=' ', names=column_names)
    sampled_stats = pa.read_csv('fess/stats/stem_stem_orientations_sampled.csv', header=None, sep=' ', names=column_names)

    max_dist = 30
    max_lat_dist = 13

    real_stats = real_stats[real_stats["dist"] < max_dist]
    real_angles = real_stats[real_stats["lateral_offset"] < max_lat_dist]["in_plane"]

    orig_real_angles = real_angles
    real_kde1 = ss.gaussian_kde(real_angles)
    real_angles = np.concatenate([-real_angles, real_angles, 2 * m.pi - real_angles])

    sampled_stats = sampled_stats[sampled_stats["dist"] < max_dist]
    sampled_angles = sampled_stats[sampled_stats["lateral_offset"] < max_lat_dist]["in_plane"]
    orig_sampled_angles = sampled_angles
    sampled_kde1 = ss.gaussian_kde(sampled_angles)
    sampled_angles = np.concatenate([-sampled_angles, sampled_angles, 2*m.pi-sampled_angles])

    real_kde = ss.gaussian_kde(real_angles, bw_method=real_kde1.factor / 4.)
    sampled_kde = ss.gaussian_kde(sampled_angles, bw_method=sampled_kde1.factor / 4.)


    print len(real_angles), len(sampled_angles)

    num_points = 200
    x_axis = np.linspace(0, 3.14, num_points)

    print sum(real_kde(x_axis)), sum(sampled_kde(x_axis))

    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    #ax.plot(x_axis, real_kde(x_axis))

    width=m.pi / (num_points * 2.)
    #ax.bar(x_axis , real_kde(x_axis), width=width, color='g')
    #ax.bar(x_axis - width, sampled_kde(x_axis), width=width, color='r')
    ax.plot(x_axis , real_kde(x_axis), color='g')
    ax.plot(x_axis , sampled_kde(x_axis), color='r')
    #ax.hist(angles, normed=True)
    #ax.hist([orig_real_angles, orig_sampled_angles], normed=True)
    #ax.hist([real_angles, sampled_angles], normed=True)

    #plt.polar(x_axis, real_kde(x_axis))
    plt.show()

if __name__ == '__main__':
    main()

