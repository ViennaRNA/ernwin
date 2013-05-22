#!/usr/bin/python

import matplotlib
import matplotlib.pyplot as plt
import pandas as pa
import math as m
import sys, numpy as np
import scipy.stats as ss
from optparse import OptionParser

import corgy.builder.energy as cbe

def main():
    usage = """
    ./python_plot_longrange.py

    Plot the longrange interactions.
    """
    num_args= 0
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    lle = cbe.LoopLoopEnergy()

    (real_data, real_d_given_i, real_d) = lle.load_data('fess/stats/temp.longrange.stats')
    (sampled_data, sampled_d_given_i, sampled_d) = lle.load_data('fess/stats/temp.longrange.stats.sampled')

    max_dist = 150
    x = np.linspace(0, 50, 100)

    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])

    ax.plot(x, real_d_given_i(x) / real_d(x), 'g', label="real")
    ax.plot(x, sampled_d_given_i(x) / sampled_d(x), 'r', label="sampled")

    ax.set_xlabel("distance")
    ax.set_ylabel("density")

    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels)
    ax.xaxis.label.set_fontsize(20)
    ax.yaxis.label.set_fontsize(20)

    for tick in ax.get_xticklabels() + ax.get_yticklabels():
        tick.set_fontsize(15)

    plt.savefig("longrange_distances.png", bbox_inches='tight')

    plt.show()

if __name__ == '__main__':
    main()

