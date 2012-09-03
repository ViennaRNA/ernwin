#!/usr/bin/python

from corgy.utilities.vector import rotation_matrix

import numpy.random as nr
import random as rand

for i in xrange(10000):
    axis = nr.random(3)
    theta = rand.random()

    r = rotation_matrix(axis, theta)
