import numpy as np
cimport numpy as np

cdef extern from "rot_mat_c.h":
    void rotation_matrix_c(char *axis, double theta, char *mat)
    void get_closer_rotation_matrix_c(char *TH, char *point, char *M, char *F, char *out_rot_mat)

ctypedef np.double_t DTYPE_t
#print dir(np)
#from math import sin, cos, sqrt

from libc.math cimport sin,cos, sqrt, acos, atan2

def dot(np.ndarray[DTYPE_t, ndim=1] v1, np.ndarray[DTYPE_t, ndim=1] v2):
    return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]

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
    rotation_matrix_c(axis.data, theta, mat.data)

    return mat

def get_closer_rotation_matrix_cython(np.ndarray[DTYPE_t, ndim=1] TH,np.ndarray[DTYPE_t, ndim=1] point,
        np.ndarray[DTYPE_t, ndim=2] M,np.ndarray[DTYPE_t, ndim=2] F,np.ndarray[DTYPE_t, ndim=2] out_rot_mat):
    get_closer_rotation_matrix_c(TH.data, point.data, M.data, F.data, out_rot_mat.data)



