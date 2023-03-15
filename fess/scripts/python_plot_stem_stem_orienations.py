#!/usr/bin/python

from __future__ import absolute_import
from __future__ import print_function
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
    parser.add_option('-s', '--sampled_stats', dest='sampled_stats', default='fess/stats/stem_stem_orientations_sampled.csv', help='The location of the sampled statistics file.')

    (options, args) = parser.parse_args()

    if len(args) < num_required_args:
        parser.print_help()
        sys.exit(1)

    column_names = ['dist', 'in_plane', 'out_of_plane', 'vec_angle', 'lateral_offset', 'ortho_offset']

    real_stats = pa.read_csv(options.real_stats, header=None, sep=' ', names=column_names)
    sampled_stats = pa.read_csv(options.sampled_stats, header=None, sep=' ', names=column_names)

    min_dist = 15 
    max_dist = 35
    max_lat_dist = 10

    real_stats = real_stats[np.logical_and(real_stats["dist"] < max_dist, real_stats["dist"] > min_dist)]

    real_angles = real_stats[real_stats["lateral_offset"] < max_lat_dist]["in_plane"]
    real_out_angles = real_stats[real_stats["lateral_offset"] < max_lat_dist]["out_of_plane"]
    real_dist = real_stats[real_stats["lateral_offset"] < max_lat_dist]["dist"]
    real_vec_angles = real_stats[real_stats["lateral_offset"] < max_lat_dist]["vec_angle"]

    orig_real_angles = real_angles
    orig_real_out_angles = real_out_angles
    real_kde1 = ss.gaussian_kde(real_angles)
    real_angles = np.concatenate([-real_angles, real_angles, 2 * m.pi - real_angles])
    real_out_angles = np.concatenate([-real_out_angles, real_out_angles, 2*m.pi - real_out_angles])

    sampled_stats = sampled_stats[np.logical_and(sampled_stats["dist"] < max_dist, sampled_stats["dist"] > min_dist)]
    sampled_angles = sampled_stats[sampled_stats["lateral_offset"] < max_lat_dist]["in_plane"]
    sampled_out_angles = sampled_stats[sampled_stats["lateral_offset"] < max_lat_dist]["out_of_plane"]
    sampled_dist = sampled_stats[sampled_stats["lateral_offset"] < max_lat_dist]["dist"]
    sampled_vec_angles = sampled_stats[sampled_stats["lateral_offset"] < max_lat_dist]["vec_angle"]

    orig_sampled_angles = sampled_angles
    orig_sampled_out_angles = sampled_out_angles
    sampled_kde1 = ss.gaussian_kde(sampled_angles)
    sampled_angles = np.concatenate([-sampled_angles, sampled_angles, 2*m.pi-sampled_angles])
    sampled_out_angles = np.concatenate([-sampled_out_angles, sampled_out_angles, 2*m.pi-sampled_out_angles])

    real_kde = ss.gaussian_kde(real_angles, bw_method=real_kde1.factor / 4.)
    real_out_kde = ss.gaussian_kde(real_out_angles, bw_method=real_kde1.factor / 4.)

    sampled_kde = ss.gaussian_kde(sampled_angles, bw_method=sampled_kde1.factor / 4.)
    sampled_out_kde = ss.gaussian_kde(sampled_out_angles, bw_method=sampled_kde1.factor / 4.)

    print(len(real_angles), len(sampled_angles))

    num_points = 200
    x_axis = np.linspace(0, 3.14, num_points)

    print(sum(real_kde(x_axis)), sum(sampled_kde(x_axis)))

    #ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    dims = (4,2)
    ax = plt.subplot2grid(dims, (0,0))
    ax1 = plt.subplot2grid(dims, (0,1))
    ax2 = plt.subplot2grid(dims, (1,0), colspan=2)
    ax3 = plt.subplot2grid(dims, (2,0), colspan=2)
    ax4 = plt.subplot2grid(dims, (3,0), colspan=1)
    ax5 = plt.subplot2grid(dims, (3,1), colspan=1)
    #ax.plot(x_axis, real_kde(x_axis))

    width=m.pi / (num_points * 2.)

    #ax.bar(x_axis , real_kde(x_axis), width=width, color='g')
    #ax.bar(x_axis - width, sampled_kde(x_axis), width=width, color='r')
    p1 = ax.plot(x_axis , real_kde(x_axis), color='g', label="real")
    p2 = ax.plot(x_axis , sampled_kde(x_axis), color='r', label="sampled")
    p3 = ax.plot(x_axis, real_kde(x_axis) -sampled_kde(x_axis), color='y', label='difference')

    p1 = ax1.plot(x_axis , real_out_kde(x_axis), color='g', label="real")
    p2 = ax1.plot(x_axis , sampled_out_kde(x_axis), color='r', label="sampled")
    p3 = ax1.plot(x_axis, real_out_kde(x_axis) -sampled_out_kde(x_axis), color='y', label='difference')
    #ax1.hist(orig_real_out_angles)
    ax2.plot(real_dist, orig_real_angles, 'go', alpha=0.5)
    ax2.plot(sampled_dist, orig_sampled_angles, 'ro', alpha=0.5)

    ax3.plot(real_dist, orig_real_out_angles, 'go', alpha=0.3)
    ax3.plot(sampled_dist, orig_sampled_out_angles, 'ro', alpha=0.3)
    #ax.hist(angles, normed=True)
    #ax.hist([orig_real_angles, orig_sampled_angles], normed=True)
    #ax.hist([real_angles, sampled_angles], normed=True)

    ax4.plot(orig_real_angles, orig_real_out_angles, 'go', alpha=0.3)
    ax4.plot(orig_sampled_angles, orig_sampled_out_angles, 'ro', alpha=0.3)

    ax5.plot(orig_real_angles, real_vec_angles, 'o', alpha=0.3)
    ax5.plot(orig_sampled_angles, sampled_vec_angles, 'o', alpha=0.3)

    ax.set_xlabel("angle")
    ax.set_ylabel("density")

    handles, labels = ax.get_legend_handles_labels()
    #ax.legend(handles, labels)
    ax.xaxis.label.set_fontsize(20)
    ax.yaxis.label.set_fontsize(20)

    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontsize(15)

    plt.savefig("stem_stem_orientations.png", bbox_inches='tight')
    #plt.polar(x_axis, real_kde(x_axis))
    plt.show()

if __name__ == '__main__':
    main()

