#!/usr/bin/python

from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import timeit, sys

from math import cos, sin, sqrt
import numpy.random as nr

from scipy import weave
from cytvec import rotation_matrix_cython

def rotation_matrix_weave(axis, theta, mat = None):
    '''
    Calculate the rotation matrix for a rotation of theta degrees around axis.

    Thanks to unutbu on StackOverflow 

    http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector

    @param axis: The axis around which to rotate
    @param theta: The angle of rotation
    @return: A matrix which can be used to perform the given rotation. The coordinates
             need only be multiplied by the matrix.
    '''
    if mat == None:
        mat = np.eye(3,3)

    support = "#include <math.h>"
    code = """
        double x = sqrt(axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]);
        double a = cos(theta / 2.0);
        double b = -(axis[0] / x) * sin(theta / 2.0);
        double c = -(axis[1] / x) * sin(theta / 2.0);
        double d = -(axis[2] / x) * sin(theta / 2.0);

        mat[0] = a*a + b*b - c*c - d*d;
        mat[1] = 2 * (b*c - a*d);
        mat[2] = 2 * (b*d + a*c);

        mat[3*1 + 0] = 2*(b*c+a*d);
        mat[3*1 + 1] = a*a+c*c-b*b-d*d;
        mat[3*1 + 2] = 2*(c*d-a*b);

        mat[3*2 + 0] = 2*(b*d-a*c);
        mat[3*2 + 1] = 2*(c*d+a*b);
        mat[3*2 + 2] = a*a+d*d-b*b-c*c;
    """

    weave.inline(code, ['axis', 'theta', 'mat'], support_code = support, libraries = ['m'])

    return mat

def rotation_matrix_numpy(axis, theta):
    mat = np.eye(3,3)
    axis = axis/sqrt(np.dot(axis, axis))
    a = cos(theta/2.)
    b, c, d = -axis*sin(theta/2.)

    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                  [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                  [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])


if __name__ == '__main__':

    setup = """
import numpy as np
import numpy.random as nr

from rotation_matrix_test import rotation_matrix_weave
from rotation_matrix_test import rotation_matrix_numpy
from rot_mat_cyt import rotation_matrix_cython

mat1 = np.eye(3,3)
theta = nr.random()
axis = nr.random(3)
"""

    print(timeit.repeat("rotation_matrix_cython(axis, theta, mat1)", setup=setup, number=100000))
    print(timeit.repeat("rotation_matrix_weave(axis, theta, mat1)", setup=setup, number=100000))
    sys.exit(1)

    print(timeit.repeat("rotation_matrix_numpy(axis, theta)", setup=setup, number=100000))
