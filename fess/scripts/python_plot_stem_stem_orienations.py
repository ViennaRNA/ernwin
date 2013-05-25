#!/usr/bin/python

import matplotlib
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
    parser.add_option('-r', '--real_stats', dest='real_stats', default='fess/stats/stem_stem_orientations.csv', help='The location of the real statistics file.')
    parser.add_option('-s', '--sampled_stats', dest='sampled_stats', default='fess/stats/stem_stem_orientations.csv', help='The location of the sampled statistics file.')

    (options, args) = parser.parse_args()

    if len(args) < num_required_args:
        parser.print_help()
        sys.exit(1)

    column_names = ['dist', 'out_of_plane', 'in_plane', 'vec_angle', 'lateral_offset', 'ortho_offset']

    real_stats = pa.read_csv(options.real_stats, header=None, sep=' ', names=column_names)
    sampled_stats = pa.read_csv(options.fake_stats, header=None, sep=' ', names=column_names)

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
    p1 = ax.plot(x_axis , real_kde(x_axis), color='g', label="real")
    p2 = ax.plot(x_axis , sampled_kde(x_axis), color='r', label="sampled")
    p3 = ax.plot(x_axis, real_kde(x_axis) -sampled_kde(x_axis), color='y', label='difference')
    #ax.hist(angles, normed=True)
    #ax.hist([orig_real_angles, orig_sampled_angles], normed=True)
    #ax.hist([real_angles, sampled_angles], normed=True)

    ax.set_xlabel("angle")
    ax.set_ylabel("density")

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)
    ax.xaxis.label.set_fontsize(20)
    ax.yaxis.label.set_fontsize(20)

    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontsize(15)

    plt.savefig("stem_stem_orientations.png", bbox_inches='tight')
    #plt.polar(x_axis, real_kde(x_axis))
    plt.show()

if __name__ == '__main__':
    main()

