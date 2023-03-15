#!/usr/bin/python

from __future__ import absolute_import
from __future__ import print_function
import sys
from optparse import OptionParser

from forgi.threedee.utilities.vector import normalize, rotation_matrix_weave
from forgi.threedee.utilities.vector import vec_angle, magnitude

from numpy import array, dot, cross, allclose, eye, ones
from numpy.random import random

import numpy as np

from math import pi, sqrt, atan2
from scipy import weave
from six.moves import range

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

def get_closer_rotation_matrix(TH, point, M, F, out_rot_mat = None):
    support = "#include <math.h>"

    code = """
    int i, j;

    double sum1 = 0;
    double TH_norm[3];
    double d1[3];

    double R[9], O[9];

    //TH_norm = normalize(TH)
    for (i = 0; i < 3; i++) 
        sum1 += TH[i] * TH[i];
    sum1 = sqrt(sum1);

    //TH_arr = array([TH_norm, TH_norm, TH_norm])
    for (i = 0; i < 3; i++)
        TH_norm[i] = TH[i] / sum1;

    for (i = 0; i < 3; i++) {
        sum1 = 0;

        for (j = 0; j < 3; j++) {
            sum1 += (M[i * 3 + j] - point[j]) * TH_norm[j];
        }
        
        d1[i] = sum1;
    }

    for ( i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            O[i*3 + j] = TH_norm[j] * d1[i] + point[j];
            R[i*3 + j] = M[i*3 + j] - O[i*3 + j];

        }
    }

    double fs_sum = 0;
    double fr_sum = 0;

    double temp_fs_sum = 0;
    double temp_fr_sum = 0;

    double n3 = 0, d3 = 0, a, r_row_sums[3], s_row_sums[3];
    double S[9], F1[9];

    for (i = 0; i < 3; i++) {
        // S = cross(R, TH)
        S[3*i + 0] = R[3*i + 1] * TH[2] - R[3*i + 2] * TH[1];
        S[3*i + 1] = R[3*i + 2] * TH[0] - R[3*i + 0] * TH[2];
        S[3*i + 2] = R[3*i + 0] * TH[1] - R[3*i + 1] * TH[0];

        s_row_sums[i] = 0;
        r_row_sums[i] = 0;

        //s_row_sums = np.sum(np.abs(S) ** 2, axis=-1) ** (1./2.)
        //r_row_sums = np.sum(np.abs(R) ** 2, axis=-1) ** (1./2.)
        //F = F - O

        temp_fs_sum = 0;
        temp_fr_sum = 0;

        for (j = 0; j < 3; j++)  {
            s_row_sums[i] += S[3*i + j] * S[3*i + j];
            r_row_sums[i] += R[3*i + j] * R[3*i + j];

            F1[3*i + j] = F[3*i + j] - O[3*i + j];

        }

        s_row_sums[i] = sqrt(s_row_sums[i]);
        r_row_sums[i] = sqrt(r_row_sums[i]);

        //np.sum(F * S, axis=-1)
        //np.sum(F * R, axis=-1)
        
        for (j = 0; j < 3; j++) {
            S[3*i + j] /= s_row_sums[i];
            R[3*i + j] /= r_row_sums[i];

            temp_fs_sum += F1[3*i + j] * S[3*i + j];
            temp_fr_sum += F1[3*i + j] * R[3*i + j];
        }

        n3 += r_row_sums[i] * temp_fs_sum;
        d3 += r_row_sums[i] * temp_fr_sum;
    }

    a = atan2(n3, d3);
    return_val = a;
    """

    support = "#include <math.h>"
    a = weave.inline(code, ['M', 'TH', 'F', 'point'], support_code = support, libraries = ['m'])

    return rotation_matrix_weave(TH, a, out_rot_mat)

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
        for i in range(1, len(moving) - 3):
            TH = (moving[i] - moving[i-1])

            rot_mat = get_closer_rotation_matrix(TH, array(moving[i-1]), array(moving[-3:]), fixed, rot_mat)

            '''
            rem_moving = array(moving[i+1:])
            rem_moving -= moving[i]
            #rem_moving = dot(rot_mat, rem_moving.transpose()).transpose()
            rem_moving = dot(rem_moving, rot_mat.transpose())
            
            rem_moving += moving[i]
            moving = moving[:i+1] + list(rem_moving)
            '''

            #distances2 = [magnitude(moving[j] - moving[j-1]) for j in range(1, len(moving))]
            #print "distances2:", distances2

            #if i == 3:
            #    sys.exit(1)

            rem_moving = moving[i+1:]
            rem_tmoving = tmoving[i+1:]

            rem_moving -= moving[i]
            #print "counter:", counter
            #print "rem_tmoving:", rem_tmoving
            dot(rem_moving, rot_mat.transpose(), rem_tmoving)
            #print "rem_moving:", rem_tmoving
            rem_tmoving += moving[i]

            tt_moving = moving
            moving = tmoving
            moving[i] = tt_moving[i]
            tmoving = tt_moving
            
            '''
            sys.exit(1)
            #dot(rot_mat, rem_moving.transpose(), rem_moving.transpose())
            rem_moving += moving[i]
            print "rem_moving:", rem_moving
            '''


        if print_d:
            rmsd = calc_rmsd(moving[-3:], fixed)
            print("iteration:", k, "rmsd:", calc_rmsd(moving[-3:], fixed)) 

    print("counter:", counter)
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

    if len(sys.argv) < 2:
        moving = ccd(moving, fixed, 10, True)
    else:
        moving = ccd(moving, fixed, iterations = int(sys.argv[1]), print_d = False)

    #print "moving:", moving

    angles2 = [vec_angle(moving[i-1] - moving[i-2], moving[i] - moving[i-1]) for i in range(2, len(moving))]
    distances2 = [magnitude(moving[i] - moving[i-1]) for i in range(1, len(moving))]

    assert(allclose(distances1, distances2))
    assert(allclose(angles1, angles2))

    #print "angles1:", angles1
    #print "angles2:", angles2


if __name__ == '__main__':
    main()

