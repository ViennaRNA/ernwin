#!/usr/bin/python

from numpy import array, dot, pi, cos, sin, cross, matrix
from numpy.linalg import inv
from numpy.testing import assert_allclose
from math import sqrt, acos, atan2
from random import random, uniform

from math import isnan

from sys import stderr

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
    return array([mult * uniform(-1, 1), mult * uniform(-1, 1), mult * uniform(-1,1)])

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


    vec1 = normalize(array(vec1))
    vec2 = normalize(array(vec2))

    if vec3 == None:
        vec3 = cross(vec1, vec2)

    vec3 = normalize(vec3)

    return array([vec1, vec2, vec3])

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
    (r,u,v) = vec

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

    @param coords: The coordinates to transform (array of n elements).
    @param new_basis: The new basis vectors (n x n matrix)
    @param old_basis: The old basis for the coordinates(n x n matrix)
    @return: The new coordinates according to the new basis
    '''
    #assert(len(coords) == len(new_basis))
    #assert(len(new_basis) == len(old_basis))

    dim = len(coords)
    #print >>stderr, "coords:", coords
    #standard_basis = get_standard_basis(dim)

    # the transform of the old basis to the standard is simply
    # equal to multiplying by the old basis

    # the transform from the standard to the new basis is
    # done by multiplying the coordinates in the standard
    # basis to the coordinates in the new_basis

    # A nice tutorial about this is located here:
    # http://tutorial.math.lamar.edu/Classes/LinAlg/ChangeOfBasis.aspx
    #print >>stderr, "old_basis:", old_basis
    standard_coords = dot(old_basis.transpose(), coords)
    #print >>stderr, "standard_coords:", standard_coords
    standard_to_new = inv(new_basis.transpose())
    #print >>stderr, "standard_to_new:", standard_to_new
    new_coords = dot(standard_to_new, standard_coords)

    return new_coords

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


def rotation_matrix(axis,theta):
    '''
    Calculate the rotation matrix for a rotation of theta degress around axis.

    Thanks to unutbu on StackOverflow 

    http://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector

    @param axis: The axis around which to rotate
    @param theta: The angle of rotation
    @return: A matrix which can be used to perform the given rotation. The coordinates
             need only be multiplied by the matrix.
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


