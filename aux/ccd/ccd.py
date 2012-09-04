#!/usr/bin/python

import sys
from optparse import OptionParser

from corgy.utilities.vector import normalize, rotation_matrix_weave
from corgy.utilities.vector import vec_angle, magnitude

from numpy import array, dot, cross, allclose, eye, ones
from numpy.random import random

import numpy as np

from math import pi, sqrt, atan2

def make_random_chain(n=12):
    """
    Return a list of random vectors, each with distance
    3.8 from the previous vector.
    """
    v=array([0,0,0])
    l=[v]
    for i in range(0, n-1):
        nv=normalize(random(3))
        nv=v + 3.8 * nv
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
    m=rotation_matrix_weave(rotation_axis, angle)
    t=-0.5*random(3)
    for i in range(-3, 0):
        v=chain[i]-t
        v=dot(v, m)
        l.append(v)
    return l

def calc_rmsd(a, b):
    assert(len(a) == len(b))

    return sqrt(sum([magnitude(a[i] - b[i]) for i in range(len(a))]) / len(a))

def point_on_line(P, A, D):
    '''
    Point on axis
    http://stackoverflow.com/questions/5227373/minimal-perpendicular-vector-between-a-point-and-a-line
    
    @param P: A point anywhere
    @param A: A point on the axis
    @param D: The direction of the axis
    '''

    D = normalize(D)
    X = A + dot((P - A), D) * D

    #print "dot((P-A), D) * D", dot((P-A), D) * D

    #print "X:", X, dot(X-P, X-A)

    return X

def get_closer_rotation_matrix(axis, point, M, F):
    TH = axis
    TH_norm = normalize(TH)
    TH_arr = array([TH_norm, TH_norm, TH_norm])
    O = point + (dot((M - point), TH_norm) * TH_arr.transpose()).transpose()
    R = M - O

    S = cross(R, TH)
    s_row_sums = np.sum(np.abs(S) ** 2, axis=-1) ** (1. / 2.)

    S = S / s_row_sums[:, np.newaxis]
    r_row_sums = np.sum(np.abs(R) ** 2, axis=-1) ** (1./ 2.) 

    R = R / r_row_sums[:, np.newaxis]
    F = array(F)
    F = F-O

    n3 = dot(np.sum(F * S, axis=-1), r_row_sums)
    d3 = dot(np.sum(F * R, axis=-1), r_row_sums)

    a = atan2(n3, d3)

    return rotation_matrix_weave(TH, a)

def ccd(moving, fixed):
    '''
    Do cyclic-coordinate descent to align the last three coordinates of moving
    to the three coordinates of fixed.
    '''
    assert(len(fixed) == 3)
    iterations = 500

    for k in xrange(iterations):
        for i in xrange(1, len(moving) - 3):
            TH = (moving[i] - moving[i-1])

            rot_mat = get_closer_rotation_matrix(TH, moving[i-1], moving[-3:], fixed[-3:])

            '''
            for j in range(i+1, len(moving)):
                moving[j] -= moving[i]
                moving[j] = dot(rot_mat, moving[j])
                moving[j] += moving[i]

            '''
            rem_moving = array(moving[i+1:])
            rem_moving -= moving[i]
            rem_moving = dot(rot_mat, rem_moving.transpose()).transpose()
            rem_moving += moving[i]
            moving = moving[:i+1] + list(rem_moving)

        rmsd = calc_rmsd(moving[-3:], fixed)
        '''
        if rmsd < 0.08:
            break
        '''

        #print "iteration:", k, "rmsd:", calc_rmsd(moving[-3:], fixed) 
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

    moving = ccd(moving, fixed)

    #print "moving:", moving

    angles2 = [vec_angle(moving[i-1] - moving[i-2], moving[i] - moving[i-1]) for i in range(2, len(moving))]
    distances2 = [magnitude(moving[i] - moving[i-1]) for i in range(1, len(moving))]

    assert(allclose(angles1, angles2))
    assert(allclose(distances1, distances2))

    #print "angles1:", angles1
    #print "angles2:", angles2


if __name__ == '__main__':
    main()

