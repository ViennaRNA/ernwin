#!/usr/bin/python

#import timeit, sys

import corgy.utilities.debug as cud

from numpy import array, dot, pi, cross, eye
from math import cos, sin
from numpy.linalg import inv, solve

import numpy as np
import math

from numpy.testing import assert_allclose
from math import acos, atan2, sqrt
from random import uniform

from math import isnan

from scipy import weave

null_array = array([0., 0., 0.])

x_array = array([1., 0., 0.])
y_array = array([0., 1., 0.])
z_array = array([0., 0., 1.])

# identity matrix
identity_matrix = array([x_array, y_array, z_array])

standard_basis = array([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
tau = 2 * pi

def get_inter_distances(vecs):
    '''
    Calculate all of the distances between the points of vecs.
    '''
    distances = []
    for i in range(len(vecs)):
        for j in range(i+1, len(vecs)):
            distances += [vec_distance(vecs[i], vecs[j])]

    return distances

def get_random_vector(mult=1.):
    return array([mult * uniform(-1, 1), mult * uniform(-1, 1), mult * uniform(-1, 1)])

def get_random_vector_pair(angle=uniform(0, pi)):
    vec1 = get_random_vector()
    vec2 = get_non_colinear_unit_vector(vec1)
    rot_vec = cross(vec1, vec2)
    rotmat = rotation_matrix(rot_vec, angle)
    vec2 = dot(rotmat, vec1)
    return (vec1, vec2)

def get_double_alignment_matrix(vp1, vp2):
    '''
    Align two sets of two vectors onto each other.

    @param vp1: A pair of two vectors.
    @param vp2: Another pair of two vectors.
    '''
    angle1 = vec_angle(vp1[0], vp1[1])
    angle2 = vec_angle(vp2[0], vp2[1])

    assert_allclose(angle1, angle2, rtol=1e-7, atol=1e-7)

    # Align the first two segments
    mat1 = get_alignment_matrix(vp1[0], vp2[0])

    # See where the second segment of the second set ends up
    # after the first alignment
    new_vp2_1 = dot(mat1, vp2[1])

    comp1 = cross(vp1[1], vp1[0])
    #comp1 = cross(vp1[0], vp1[1])
    comp2 = cross(vp1[0], comp1) # should be along the plane of vp1[0] and vp1[1]

    basis1 = create_orthonormal_basis(normalize(vp1[0]), normalize(comp2))
    rej2 = change_basis(new_vp2_1, basis1, standard_basis)

    angle = atan2(rej2[2], rej2[1])

    mat2 = rotation_matrix(vp1[0], angle)

    #return dot(mat1, mat2)
    return dot(mat2, mat1)

def get_alignment_matrix(vec1, vec2):
    '''
    Return a rotation matrix that will align vec1 along vec2.

    @param vec1: The target vector.
    @param vec2: The vector to be aligned.
    '''

    comp = cross(vec1, vec2)
    angle = vec_angle(vec1, vec2)

    return rotation_matrix(comp, angle)


def get_non_colinear_unit_vector(vec): 
    '''
    Get a unit vector that does not lie on the line defined by vec.

    This is done by creating a vector along the least represented axis in vec.

    @param vec: The vector under consideration.
    @return: A vector along an axis.
    '''
    absvec = [abs(v) for v in vec]
    m = min(absvec)
    ind = absvec.index(m) 
    unit = [0., 0., 0.] 
    unit[ind] = 1. 
           
    return array(unit)

def create_orthonormal_basis1(vec1, vec2=None, vec3=None):
    '''
    Create an orthonormal basis using the provided vectors.

    If more than one is provided, it must be orthogonal to 
    the others.
    '''
    if vec2 == None:
        vec2 = get_non_colinear_unit_vector(vec1)
    #else:
    #    assert_allclose(dot(vec2, vec1), 0., rtol=1e-7, atol=1e-7)


    vec1 = normalize(array(vec1))
    vec2 = normalize(array(vec2))

    if vec3 == None:
        vec3 = cross(vec1, vec2)

    vec3 = normalize(vec3)

    return array([vec1, vec2, vec3])


def create_orthonormal_basis(vec1, vec2=None, vec3=None):
    '''
    Create an orthonormal basis using the provided vectors.

    If more than one is provided, it must be orthogonal to 
    the others.
    '''
    if vec2 == None:
        vec2 = get_non_colinear_unit_vector(vec1)
    #else:
    #    assert_allclose(dot(vec2, vec1), 0., rtol=1e-7, atol=1e-7)

    if vec3 == None:
        vec3 = cross(vec1, vec2)

    vec1 /= magnitude(vec1)
    vec2 /= magnitude(vec2)
    vec3 /= magnitude(vec3)

    return array([vec1, vec2, vec3])

def time_cob1():
    vec1 = get_random_vector()
    vec2 = get_random_vector()

    basis = create_orthonormal_basis1(vec1, vec2)

def time_cob2():
    vec1 = get_random_vector()
    vec2 = get_random_vector()

    basis = create_orthonormal_basis(vec1, vec2)

def time_cob():
    t1 = timeit.Timer("time_cob1()", "from corgy.utilities.vector import time_cob1")
    t2 = timeit.Timer("time_cob2()", "from corgy.utilities.vector import time_cob2")

    print t1.repeat(number=10000)
    print t2.repeat(number=10000)

def spherical_cartesian_to_polar(vec):
    '''
    Return a parameterization of a vector of 3 coordinates:

    x = r sin u cos v
    y = r sin u sin v
    z = r cos u

    0 <= u <= pi
    -pi <= v <= pi

    Where u is the polar angle and v is the azimuth angle.

    @param vec: A vector of 3 cartesian coordinates.
    @return: (r, u, v)
    '''
    r = magnitude(vec)
    u = acos(vec[2] / r)
    v = atan2(vec[1], vec[0])

    assert_allclose(vec[0], r * sin(u) * cos(v), rtol=1e-7, atol=1e-7)
    return (r, u, v)

def spherical_polar_to_cartesian(vec):
    '''
    Convert spherical polar coordinates to cartesian coordinates:

    See the definition of spherical_cartesian_to_polar.

    @param vec: A vector of the 3 polar coordinates (r, u, v)
    @return: (x, y, z)
    '''
    (r, u, v) = vec

    x = r * sin(u) * cos(v)
    y = r * sin(u) * sin(v)
    z = r * cos(u)

    return [x, y, z]

def get_standard_basis(dim):
    '''
    Get a standard basis for the given dimension.

    For 2D, this equals [[1.,0.],[0.,1.]]

    @param dim: The dimension of the vector space.
    @return: A vector of vectors that constitute the standard basis.
    '''

    standard_basis = [[0. for j in range(dim)] for i in range(dim)]
    for i in range(dim):
        standard_basis[i][i] = 1.
    standard_basis = array(standard_basis)

    return standard_basis

def change_basis(coords, new_basis, old_basis):
    '''
    Change the basis of coordinates to a new basis. In a regular structure
    we have the coordinates in the regular cartesian coordinate system. For helix-helix
    orientations, however, we want to express the coordinates in a coordinate system
    defined by the first helix.

    The new basis will consist of the axis of the first helix, one of its twist elements
    (the one closest to the second vector) and a third vector orthogonal to the previous
    two.

    # http://tutorial.math.lamar.edu/Classes/LinAlg/ChangeOfBasis.aspx

    @param coords: The coordinates to transform (array of n elements).
    @param new_basis: The new basis vectors (n x n matrix)
    @param old_basis: The old basis for the coordinates(n x n matrix)
    @return: The new coordinates according to the new basis
    '''
    #assert(len(coords) == len(new_basis))
    #assert(len(new_basis) == len(old_basis))

    dim = len(coords)
    #print "coords:", coords
    standard_coords = dot(old_basis.transpose(), coords)
    '''
    #print "standard_coords:", standard_coords
    standard_to_new = inv(new_basis.transpose())
    #print "standard_to_new:", standard_to_new
    new_coords = dot(standard_to_new, standard_coords)
    print "new_coords:", new_coords
    '''

    new_coords = solve(new_basis.transpose(), standard_coords)
    #print "new_coords1:", new_coords1

    return new_coords


def change_basis1(coords, new_basis, old_basis):
    '''
    '''

    dim = len(coords)
    standard_coords = dot(old_basis.transpose(), coords)
    standard_to_new = inv(new_basis.transpose())
    new_coords = dot(standard_to_new, standard_coords)

    return new_coords


def change_basis2(coords, new_basis, old_basis):
    '''
    '''

    dim = len(coords)
    standard_coords = dot(old_basis.T, coords)
    new_coords = solve(new_basis.T, standard_coords)

    return new_coords

def change_basis1_benchmark():
    coords = get_random_vector(10.)
    basis1 = array([get_random_vector(10.) for i in range(3)])
    basis2 = array([get_random_vector(10.) for i in range(3)])

    nc = change_basis1(coords, basis1, basis2)

def change_basis2_benchmark():
    coords = get_random_vector(10.)
    basis1 = array([get_random_vector(10.) for i in range(3)])
    basis2 = array([get_random_vector(10.) for i in range(3)])

    nc = change_basis2(coords, basis1, basis2)

def time_basis1():
    t1 = timeit.Timer("change_basis1_benchmark()","from corgy.utilities.vector import change_basis1_benchmark")
    print t1.repeat(number=10000)

def time_basis2():
    t2 = timeit.Timer("change_basis2_benchmark()","from corgy.utilities.vector import change_basis2_benchmark")
    print t2.repeat(number=10000)

def vector_rejection(a, b):
    '''
    Return the vector rejection of a from b. In other words, return the orthogonal
    projection of a onto the plane orthogonal to b.

    @param a: The vector to be projected.
    @param b: The vector defining the normal of the plane.
    @return: The rejection of the vector a from b. (a - (dot(a, b) / dot(b, b)) * b)
    '''

    n = dot(a, b)
    d = dot(b, b)
    return a - (n / d) * b

def rotation_matrix_weave(axis, theta, mat = None):
    '''
    Calculate the rotation matrix for a rotation of theta degrees around axis.
    
    Implemented in C++ using the weave module. Runs approximately 6x faster than 
    the numpy version below if no mat is passed in and around 20x faster if mat is
    passed in.

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


def vector_set_rmsd(set1, set2):
    '''
    Calculate the rmsd between two sets of vectors.

    @param set1: A matrix
    @param set2: Another matrix.

    @return: The rmsd between the rows of the matrix.
    '''
    rmsd = 0
    count = 0
    for i in range(len(set1)):
        rmsd += magnitude(set2[i] - set1[i])
        count += 1
    rmsd /= count
    return sqrt(rmsd)

def rotation_matrix(axis, theta):
    '''
    Calculate the rotation matrix for a rotation of theta degrees around axis.

    Thanks to unutbu on StackOverflow 

    http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector

    @param axis: The axis around which to rotate
    @param theta: The angle of rotation
    @return: A matrix which can be used to perform the given rotation. The coordinates
             need only be multiplied by the matrix.
    '''
    axis = axis/sqrt(dot(axis, axis))
    a = cos(theta/2)
    b, c, d = -axis*sin(theta/2)

    return array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                  [2*(b*c+a*d), a*a+c*c-b*b-d*d, 2*(c*d-a*b)],
                  [2*(b*d-a*c), 2*(c*d+a*b), a*a+d*d-b*b-c*c]])

def get_vector_centroid(crds1):
    '''
    Find the centroid of a set of vectors.

    @param crds: A matrix containing all of the vectors.

    @return: The centroid of the rows of the matrix crds.
    '''
    centroid1 = array([0., 0., 0.])

    for i in range(len(crds1)):
        centroid1 += crds1[i]

    centroid1 /= float(len(crds1))

    for i in centroid1:
        if isnan(i):
            raise Exception('nan encountered')

    return centroid1

def center_on_centroid(crds1):
    centroid1 = get_vector_centroid(crds1)
    crds = []

    for i in range(len(crds1)):
        crds += [crds1[i] - centroid1]

    return array(crds)

def magnitude(vec):
    '''
    Return the magnitude of a vector (|V|).

    @param vec: The vector in question.
    @return: The magnitude of the vector.
    '''

    return sqrt(dot(vec, vec))

def time_mag1():
    vec1 = get_random_vector()

    return sqrt(dot(vec1, vec1))

def time_mag2():
    vec1 = get_random_vector()

    return sqrt(dot(vec1, vec1))

def time_mag():
    t1 = timeit.Timer("time_mag1()", "from corgy.utilities.vector import time_mag1")
    t2 = timeit.Timer("time_mag2()", "from corgy.utilities.vector import time_mag2")

    print t1.repeat(number=10000)
    print t2.repeat(number=10000)


def normalize(vec):
    '''
    Normalize a vector so that its magnitude becomes 1.0 while
    its direction remains the same.

    @param vec: The vector in question.
    @return: A normalized version of the vector.
    '''

    return vec / magnitude(vec)

def vec_angle(vec1, vec2):
    '''
    Get the angle between two vectors using the identity:

    A * B = |A||B| cos t

    Where A and B are two vectors and t is the angle between them.

    @param vec1: The first vector (A)
    @param vec2: The second vector (B)
    @return: The angle between the two vectors.
    '''

    vec1n = normalize(vec1)
    vec2n = normalize(vec2)

    cud.pv('dot(vec1n, vec2n)')
    angle = acos(dot(vec1n, vec2n))
    return angle

def vec_dot(a, b):
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]

'''
def cross(a, b):
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]

    return c
'''

def vec_distance(vec1, vec2):
    return sqrt(dot(vec2 - vec1, vec2 - vec1))


def line_segment_distance(s1_p0, s1_p1, s2_p0, s2_p1):
    '''
    Calculate the two points on each of the segments that are closest to
    each other. The first segment is defined as p1->p2 and the second as
    p3->p4.

    Code shamelessly translated from:
    http://softsurfer.com/Archive/algorithm_0106/algorithm_0106.htm#dist3D_Segment_to_Segment

    @param s1_p0: The start of the first segment
    @param s1_p1: The end of the first segment

    @param s2_p0: The start of the second segment
    @param s2_p1: The end of the second segment
    @return: A tuple of points (i1,i2) containing the point i1 on s1
        closest to the point i2 on segment s2
    '''
    u = s1_p1 - s1_p0
    v = s2_p1 - s2_p0
    w = s1_p0 - s2_p0
    
    a = dot(u,u)        # always >= 0
    b = dot(u,v)
    c = dot(v,v)        # always >= 0
    d = dot(u,w)
    e = dot(v,w)
    D = a*c - b*b       # always >= 0
    sD = D      # sc = sN / sD, default sD = D >= 0
    tD = D      # tc = tN / tD, default tD = D >= 0

    SMALL_NUM = 0.000001

    # compute the line parameters of the two closest points
    if (D < SMALL_NUM):  # the lines are almost parallel
        sN = 0.0        # force using point P0 on segment S1
        sD = 1.0        # to prevent possible division by 0.0 later
        tN = e
        tD = c
    else:                # get the closest points on the infinite lines
        sN = (b*e - c*d)
        tN = (a*e - b*d)
        if (sN < 0.0):      # sc < 0 => the s=0 edge is visible
            sN = 0.0
            tN = e
            tD = c
        elif (sN > sD):  # sc > 1 => the s=1 edge is visible
            sN = sD
            tN = e + b
            tD = c

    if (tN < 0.0):           # tc < 0 => the t=0 edge is visible
        tN = 0.0
        # recompute sc for this edge
        if (-d < 0.0):
            sN = 0.0
        elif (-d > a):
            sN = sD
        else:
            sN = -d
            sD = a
    elif (tN > tD):      # tc > 1 => the t=1 edge is visible
        tN = tD
        # recompute sc for this edge
        if ((-d + b) < 0.0):
            sN = 0
        elif ((-d + b) > a):
            sN = sD
        else:
            sN = (-d + b)
            sD = a
            
    # finally do the division to get sc and tc
    sc = 0.0 if abs(sN) < SMALL_NUM else sN / sD
    tc = 0.0 if abs(tN) < SMALL_NUM else tN / tD

    # get the difference of the two closest points
    dP = w + (sc * u) - (tc * v)  # = S1(sc) - S2(tc)

    '''
    cud.pv('u')
    cud.pv('v')
    cud.pv('sc')
    cud.pv('tc')

    cud.pv('sc')
    cud.pv('tc')

    cud.pv('s1_p0 + sc * u')
    cud.pv('s2_p0 + tc * v')

    cud.pv('dP')
    '''

    return (s1_p0 + sc * u, s2_p0 + tc * v)
