#!/usr/bin/python

import sys
from numpy import dot, array
from shape_fit import rotation_matrix, vec_angle
from math import pi

vec1 = array([1,4,0])
vec2 = array([1,4,0])
mats = []

increment = 100

for i in range(increment):
    rot_angle = 4 * pi * float(i+1) * (1 / float(increment))
    mat = rotation_matrix(array([0., 1., 0.]), rot_angle)
    vec3 = dot(mat, vec1)
    print >>sys.stderr, "vec3:", vec3
    angle = vec_angle(vec2, vec3)
    print rot_angle, angle
