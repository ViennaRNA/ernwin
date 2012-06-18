#!/usr/bin/python

from numpy import array, dot, pi, cos, sin, cross
from math import sqrt, acos
from random import random

tau = 2 * pi

def get_non_colinear_unit_vector(vec): 
    m = min(vec) 
    ind = list(vec).index(m) 
    unit = [0., 0., 0.] 
    unit[ind] = 1. 
           
    return array(unit)

def rotation_matrix(axis,theta):
    '''
    Calculate the rotation matrix for a rotation of theta degress around axis.

    Thanks to unutbu on StackOverflow 

    http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    '''
    axis = axis/sqrt(dot(axis,axis))
    a = cos(theta/2)
    b,c,d = -axis*sin(theta/2)
    return array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                  [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                  [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])


def get_vector_centroid(crds1):
    centroid1 = array([0., 0., 0.])

    for i in range(len(crds1)):
        centroid1 += crds1[i]

    centroid1 /= float(len(crds1))

    return centroid1

def center_on_centroid(crds1):
    centroid1 = get_vector_centroid(crds1)
    crds = []

    for i in range(len(crds1)):
        crds += [crds1[i] - centroid1]

    return array(crds)

def magnitude(vec):
    return sqrt(dot(vec, vec))

def normalize(vec):
    return vec / magnitude(vec)

def vec_angle(vec1, vec2):
    vec1n = normalize(vec1)
    vec2n = normalize(vec2)

    #print >>sys.stderr, "vec1n:", vec1n
    #print >>sys.stderr, "vec2n:", vec2n
    #print >>sys.stderr, "dot(v1, v2):", dot(vec1n, vec2n)

    angle = acos(dot(vec1n, vec2n))
    return angle

def dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

def cross(a, b):
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]

    return c
def vec_distance(vec1, vec2):
    return sqrt(dot(vec2 - vec1, vec2 - vec1))


