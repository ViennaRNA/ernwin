#!/usr/bin/python

from __future__ import absolute_import
import sys
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser

def main():
    usage = """
    python plot_two_columns.py [file]

    Plot a file containing two columns, containing the x and y values,
    respectively.
    """
    num_args= 1
    parser = OptionParser(usage=usage)

    parser.add_option('-o', '--output-file', dest='output_file', default=None, help="Save the plot to a file.", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    data = np.genfromtxt(args[0], delimiter=' ')

    plt.plot(data[:,0], data[:,1])

    if options.output_file is not None:
        plt.savefig(options.output_file)

    plt.show()

if __name__ == '__main__':
    main()

