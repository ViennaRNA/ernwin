#!/usr/bin/python

from __future__ import absolute_import
from __future__ import print_function
import sys
from optparse import OptionParser

from cytvec import rotation_matrix_cython
from cytvec import normalize_in_place
from cytvec import magnitude
from cytvec import vec_angle
from cytvec import get_closer_rotation_matrix_cython
from cytvec import ccd_cython

from numpy import array, dot, cross, allclose, eye, ones
from numpy.random import random, randint

import numpy as np

from scipy import weave

from math import pi, sqrt, atan2
from six.moves import range

def make_random_chain(n=12):
    """
    Return a list of random vectors, each with distance
    3.8 from the previous vector.
    """
    v=array([0,0,0])
    l=[v]
    for i in range(0, n-1):
        nv = random(3)
        normalize_in_place(nv)
        nv=v + 38 * nv
        l.append(nv)
        v=nv

    return l

def rotate_last_three(chain):
    """
    Take the last three vectors of chain, copy them,
    and apply a random rotation and translation.
    """
    l=[]
    rotation_axis=random(3)
    angle=0.1+0.1*pi*random()
    m = np.eye(3,3)
    rotation_matrix_cython(rotation_axis, angle, m)
    t=-0.5*random(3)
    for i in range(-3, 0):
        v=chain[i]-t
        v=dot(v, m)
        l.append(v)
    return l

def calc_rmsd(a, b):
    assert(len(a) == len(b))

    return sqrt(sum([magnitude(a[i] - b[i]) for i in range(len(a))]) / len(a))

def ccd(moving, fixed, iterations=10, print_d=True):
    '''
    Do cyclic-coordinate descent to align the last three coordinates of moving
    to the three coordinates of fixed.
    '''
    assert(len(fixed) == 3)
    fixed = array(fixed)
    counter = 0
    moving1 = array(moving)
    moving2 = array(moving)
    tmoving = moving2
    moving = moving1
    rot_mat = np.eye(3,3)

    for k in range(iterations):
        prev_i = 1
        for i in range(1, len(moving) - 3, 1):
            TH = (moving[i] - moving[i-1])

            get_closer_rotation_matrix_cython(TH, array(moving[i-1]), array(moving[-3:]), fixed, rot_mat)

            #distances1 = [magnitude(moving[j] - moving[j-1]) for j in range(1, len(moving))]
            rem_moving = moving[i+1:]
            rem_tmoving = tmoving[i+1:]

            '''
            if i == 2:
                print "before:", moving[i+1]
                print "before translated:", moving[i+1] - moving[i]
            '''

            rem_moving -= moving[i]

            dot(rem_moving, rot_mat.transpose(), rem_tmoving)

            '''
            if i == 2:
                print "after pre-translated:", rem_tmoving[0]
            '''

            rem_tmoving += moving[i]
            moving[:] = tmoving

            '''
            if i == 2:
                print "rot_mat:", rot_mat
                print "moving[i+1]", moving[i+1]

            if i == 2:
                return
            '''
            #distances2 = [magnitude(moving[j] - moving[j-1]) for j in range(1, len(moving))]
            #assert(np.allclose(distances1, distances2))

            prev_i = i


        rmsd = calc_rmsd(moving[-3:], fixed)
        if print_d:
            print("iteration:", k, "rmsd:", calc_rmsd(moving[-3:], fixed)) 

    print("counter:", counter, "rmsd_d:", rmsd)
    return moving

def main():
    # Moving segment
    moving=make_random_chain(20)
    # Fixed segment 
    # Last three residues of the moving segment
    # after applying a random rotation/translation
    fixed=rotate_last_three(moving)

    angles1 = [vec_angle(moving[i-1] - moving[i-2], moving[i] - moving[i-1]) for i in range(2, len(moving))]
    distances1 = [magnitude(moving[i] - moving[i-1]) for i in range(1, len(moving))]

    #print "moving:", moving

    moving = array(moving)
    fixed = array(fixed)
    #print "moving[-3]:", array(moving)[-3:]
    points = array([i for i in range(1, len(moving) - 3, 2)])
    if len(sys.argv) < 2:
        #ccd(array(moving), array(fixed), 20)
        print("=================")
        ccd_cython(moving, fixed, points, len(moving)-3, 20)
    else:
        ccd_cython(moving, fixed, points, len(moving)-3, int(sys.argv[1]))

    print(calc_rmsd(moving[-3:], fixed))

    '''
    if len(sys.argv) < 2:
        moving = ccd(moving, fixed, 10, True)
    else:
        moving = ccd(moving, fixed, iterations = int(sys.argv[1]), print_d = False)
    '''

    #print "moving:", moving

    angles2 = [vec_angle(moving[i-1] - moving[i-2], moving[i] - moving[i-1]) for i in range(2, len(moving))]
    distances2 = [magnitude(moving[i] - moving[i-1]) for i in range(1, len(moving))]

    assert(allclose(distances1, distances2))
    assert(allclose(angles1, angles2))

    #print "angles1:", angles1
    #print "angles2:", angles2


if __name__ == '__main__':
    main()

