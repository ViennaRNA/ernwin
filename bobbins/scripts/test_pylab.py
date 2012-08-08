#!/usr/bin/python

import sys

from pylab import plot, show, draw, ion, ioff
from random import random
from time import sleep

from scipy.stats import norm, gaussian_kde

import matplotlib.pyplot as plt
import numpy as np


def test1():
    fig = plt.figure(figsize=(9,9))
    ax1 = fig.add_subplot(2,2,1)
    ax2 = fig.add_subplot(2,2,2)

    ax3 = fig.add_subplot(2,2,3)

    for i in range(20):
        if random() < 0.3:
            #pl.figure(0)
            ax1.plot(random(), random(), 'bo')
        elif random() < 0.3:
            ax3.plot(random(), random(), 'go')
        else:
            #pl.figure(1)
            ax2.plot(random(), random(), 'ro')

        draw()

def main():
    if len(sys.argv) < 0:
        print >>sys.stderr, "Usage: usage"
        sys.exit(1)

    plt.ion()

    fig = plt.figure(figsize=(9,9))
    ax = fig.add_subplot(1,1,1)

    m1 = norm.rvs(0, 1, size=20)
    m2 = norm.rvs(0, 1, size=20)

    xmin = m1.min()
    xmax = m1.max()

    ymin = m2.min()
    ymax = m2.max()

    X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([X.ravel(), Y.ravel()])
    values = np.vstack([m1, m2])
    kernel = gaussian_kde(values)
    Z = np.reshape(kernel(positions).T, X.shape)

    ax.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,
                      extent=[xmin, xmax, ymin, ymax])

    ax.plot(m1, m2, 'k.', markersize=2)
    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    plt.show()

    ioff()
    show()

if __name__ == '__main__':
    main()

