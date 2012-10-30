import numpy as np
cimport numpy as np

cdef extern from "rot_mat_c.h":
    void rotation_matrix_c(double *axis, double theta, double *mat)
    void get_closer_rotation_matrix_c(double *TH, double *point, double *M, double *F, double *out_rot_mat)
    void ccd_c(double *moving, int len_moving, double *fixed, long *points, int num_points, int moving_end, int iterations)
    double c_dot(double *v1, double *v2)

ctypedef np.double_t DTYPE_t
#print dir(np)
#from math import sin, cos, sqrt

from libc.math cimport sin,cos, sqrt, acos, atan2, pow

def dot(np.ndarray[DTYPE_t, ndim=1] v1, np.ndarray[DTYPE_t, ndim=1] v2):
    #print "yo yo"
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]

def dot(np.ndarray[DTYPE_t, ndim=1] v1, np.ndarray[DTYPE_t, ndim=1] v2):

def magnitude(np.ndarray[DTYPE_t, ndim=1] vec):
    cdef double x = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2])

    return x

def vec_angle(np.ndarray[DTYPE_t, ndim=1] v1, np.ndarray[DTYPE_t, ndim=1] v2):
    cdef double m1 = magnitude(v1)
    cdef double m2 = magnitude(v2)

    cdef double dot = v1[0] * v2[0] / (m1 * m2) + v1[1] * v2[1] / (m1 * m2) + v1[2] * v2[2] / (m1 * m2)
    cdef double angle = acos(dot)

    return angle

def normalize_in_place(np.ndarray[DTYPE_t, ndim=1] vec):
    cdef double x = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2])
    vec[0] = vec[0] / x
    vec[1] = vec[1] / x
    vec[2] = vec[2] / x

#cimport cython
#@cython.boundscheck(False) # turn of bounds-checking for entire function
def rotation_matrix_cython(np.ndarray[DTYPE_t, ndim=1] axis, double theta, np.ndarray[DTYPE_t, ndim=2] mat):
    #mat = np.eye(3,3)
    #print dir(axis)
    rotation_matrix_c(<double *> axis.data, theta, <double *>mat.data)

    return mat

def get_closer_rotation_matrix_cython(np.ndarray[DTYPE_t, ndim=1] TH,np.ndarray[DTYPE_t, ndim=1] point,
        np.ndarray[DTYPE_t, ndim=2] M, np.ndarray[DTYPE_t, ndim=2] F,np.ndarray[DTYPE_t, ndim=2] out_rot_mat):
    get_closer_rotation_matrix_c(<double *> TH.data, <double *> point.data, <double *> M.data, <double *> F.data, <double *>out_rot_mat.data)

def ccd_cython(np.ndarray[DTYPE_t, ndim=2] moving, np.ndarray[DTYPE_t, ndim=2] fixed, np.ndarray[np.int64_t, ndim=1] points, int moving_end, int iterations=10):
    cdef int len_points = len(points)
    cdef int len_moving = len(moving)
    ccd_c(<double *> moving.data, len_moving, <double *> fixed.data, <long *> points.data, len_points, moving_end, iterations) 

