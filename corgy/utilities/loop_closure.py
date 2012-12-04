#!/usr/bin/python

from corgy.utilities.data_structures import DefaultDict
from corgy.utilities.vector import vec_distance, normalize, rotation_matrix, vec_angle, magnitude
from corgy.utilities.vector import get_non_colinear_unit_vector

import sys
from random import random, normalvariate, seed
from math import cos, sin, pi, sqrt, acos, pi
from numpy import array, dot, cross

def print_segment(*karg, **kwargs):
    pass

def connect_to_start(prev_x, prev_y, prev_z, lengths):
    """
    Position the remaining two segments so that they connect
    to the origin.
    """

    r0 = lengths[0]
    r1 = lengths[1]

    x0 = 0.
    y0 = 0.

    x1 = prev_x
    y1 = prev_y

    d = sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)

    if d > r0 + r1:
        print >> sys.stderr, "Last two segments are not long enough to make the connection between the two points."
        sys.exit(1)
    if d < abs(r0 - r1):
        print >> sys.stderr, "Concentric circles."
        sys.exit(1)

    a = (r0 ** 2 - r1 ** 2 + d ** 2) / (2 * d)
    #print >>sys.stderr, "d: %f sum(lengths): %f r0: %f r1: %f a: %f" % (d, sum(lengths), r0, r1, a )
    #print >>sys.stderr, "r0 ** 2 - a ** 2: %f, r1 ** 2 - a ** 2: %f" % ( r0 ** 2 - a ** 2, r1 ** 2 - a ** 2)
    h = sqrt(r0 ** 2 - a ** 2)

    x2 = x0 + a * (x1 - x0) / d
    y2 = y0 + a * (y1 - y0) / d

    x3_1 = x2 + h * (y1 - y0) / d
    x3_2 = x2 - h * (y1 - y0) / d

    y3_1 = y2 - h * (x1 - x0) / d
    y3_2 = y2 + h * (x1 - x0) / d

    #print >>sys.stderr, "r0: %f, r1: %f, d: %f, a: %f h: %f" % (r0, r1, d, a, h)
    #print >>sys.stderr, "x0: %f, y0: %f, x1: %f, y1: %f, x2: %f, y2: %f, x3_1: %f, y3_1: %f, x3_2: %f, y3_2: %f" % (x0, y0, x1, y1, x2, y2, x3_1, y3_1, x3_2, y3_2)

    print_segment(x0, y0, 0., x3_1, y3_1, 0., color='red')
    #print_segment(x0, y0, x3_2, y3_2, color='red')

    print_segment(x1, y1, 0., x3_1, y3_1, 0., color='red')
    #print_segment(x1, y1, x3_2, y3_2, color='red')

def IntersectSpheres(o1, r1, o2, r2):
	"""
	Sourced from:
	
	http://python.rhino3d.com/threads/389-IntersectSpheres
	
	Compute the intersection of two spheres.
	Parameters:
		o1 = Center of the first sphere.
		r1 = Radius of the first sphere.
		o2 = Center of the second sphere.
		r2 = Radius of the second sphere
	Returns:	
                Array    An array of intersection results, if successful. The results are as follows:
                Element    Type         Description
                                0              Number     The type of intersection, where 0 = point, 1 = circle, 2 = spheres are identical.
                                1              Array	        If a point intersection, then an array identifying the 3-D intersection location.
                                                                If a circle intersection, then the circle's plane. The origin of the plane will be 
                                                                the center point of the circle.
                                2              Number     If a circle intersection, then the radius of the circle. (circle intersection only)
		Null If not successful, or on error
	"""
	
	# distance [o1,o2]
	d = sqrt( (o2[0] - o1[0])**2 + (o2[1] - o1[1])**2 + (o2[2] - o1[2])**2 )

	if ( d==0 ):
		if ( r1==r2 ): return [2]  # identical shperes
		else: return  # none secant spheres

	elif ( d>0 and d<abs(r2 - r1) ): return  # none secant spheres

	elif ( d==abs(r2 - r1) ):
		point = [0,0,0]
		for i in range(3):
			point[i] = o1[i] + r1 * (o2[i] - o1[i]) / d
		return [0, point]  # interior tangent spheres

	elif ( d > abs(r1 - r2) and d < abs(r1 + r2) ):
		point = [0,0,0]
		n = [0,0,0]
		x1 = abs(d**2 + r1**2 - r2**2) / (2 * d)
		x2 = abs(d**2 + r2**2 - r1**2) / (2 * d)
		r = sqrt(r1**2 - x1**2)
		for i in range(3):
			n[i] = (o2[i] - o1[i]) / d
			point[i] = o1[i] + x1 * (o2[i] - o1[i]) / d
		return [1, point, n, r]  # intersection is a circle

	elif ( d==abs(r1 + r2) ):
	    point = [0,0,0]
	    for i in range(3):
			point[i] = o1[i] + r1 * (o2[i] - o1[i]) / d
	    return [0, point]  # exterior tangent spheres
	
	elif ( d>abs(r1 + r2) ): return  # none secant spheres
	
def get_point_on_circle(o, n, r):
    #print >>sys.stderr, "origin length:", sqrt(sum([i ** 2 for i in o]))
    #print >>sys.stderr, "normal length:", sqrt(sum([i ** 2 for i in n]))
    #print >>sys.stderr, "radius:", r

    #print >>sys.stderr, "o:", o
    #print >>sys.stderr, "n:", n

    n = normalize(n)
    unit = normalize(get_non_colinear_unit_vector(n))

    U = normalize(array(cross(unit, n)))
    V = normalize(array(cross(n, U)))

    #print >>sys.stderr, "U:", U
    #print >>sys.stderr, "V:", V
    #print >>sys.stderr, "dot(U,V):", dot(U,V)
    #print >>sys.stderr, "dot(U,n):", dot(U,n)
    #print >>sys.stderr, "dot(V,n):", dot(V,n)

#    t = random() * 2 * pi
    t = 0
    p = o + (r * cos(t)) * U + (r * sin(t)) * V
    #print >>sys.stderr, "p:", p

    return p

def connect_in_3d(o1, o2, lengths):
    """
    Connect two points in 3d
    """
    r1 = lengths[0]
    r2 = lengths[1]

    #print >>sys.stderr, "r1:", r1, "r2:", r2
    #print >>sys.stderr, "o1:", o1, "o2:", o2
    #print >>sys.stderr, "dist(o1 to origin)", vec_distance(o1, array([0., 0., 0.]))

    #intersection = IntersectSpheres(o1, r1, o2, r2) 
    intersection = IntersectSpheres(o2, r2, o1, r1) 

    if intersection == None:
        print >>sys.stderr, "Spheres do not intersect."
        sys.exit(1)

    if intersection[0] == 2:
        print >>sys.stderr, "Identical spheres"
        sys.exit(1)

    if intersection[0] == 0:
        # Intersection is a point 
        print_segment(o1, intersection[1], color='red')
        print_segment(o2, intersection[1], color='red')

    elif intersection[0] == 1:
        #print >>sys.stderr, "plane:", intersection

        o = array(intersection[1])
        n = array(intersection[2])
        r = intersection[3]

        p = get_point_on_circle(o, n, r)

        print_segment(o1, p, color='red')
        print_segment(o2, p, color='red')

        print >>sys.stderr, "dist1:", vec_distance(o1, p)
        print >>sys.stderr, "dist2:", vec_distance(o2, p)


def get_random_spherical1(d):
    x = normalvariate(0, 1)
    y = normalvariate(0, 1)
    z = normalvariate(0, 1)

    p = d / sqrt(x ** 2 + y ** 2 + z ** 2)

    return (x * p, y * p, z * p)

def get_random_spherical_by_rotation(pos, dist, min_angle, max_angle):
    if pos[0] == 0. and pos[1] == 0. and pos[2] == 0.:
        pos = array([1., 0., 0.])
    
    ov = normalize(pos)
    uv1 = cross(get_non_colinear_unit_vector(pos), pos)
    uv2 = cross(pos, uv1)

    angle1 = min_angle + random() * (max_angle - min_angle)
    angle2 = random() * 2 * pi

    mat1 = rotation_matrix(uv1, angle1)
    mat2 = rotation_matrix(ov, angle2)

    nv1 = dot(mat1, ov)
    nv2 = normalize(dot(mat2, nv1))

    print >>sys.stderr, "angle1:", angle1

    return dist * nv2

def get_random_spherical(r, u = None, v = None):
    if u == None:
        u = random() * pi
    if v == None:
        v = random() * 2. * pi

    x = r * sin(u) * cos(v)
    y = r * sin(u) * sin(v)
    z = r * cos(u)

    return (x, y, z)

def get_triangle_angle(a, b, c):
    print >>sys.stderr, "a:", a
    print >>sys.stderr, "b:", b
    print >>sys.stderr, "c:", c

    try:
        return acos((a ** 2 + b ** 2 - c ** 2) / (2 * a * b))
    except ValueError as e:
        print >>sys.stderr, e
        print >>sys.stderr, "a: %f b: %f c: %f\n" % (a, b, c)
        sys.exit(1)

def get_angle_constraint(curr_pos, curr_seg, distance):
    intersection = IntersectSpheres((0., 0., 0.), distance, curr_pos, curr_seg)

    if intersection == None:
        #spheres do not intersect
        #therefore no lower limit on the angle
        min_angle = 0.
    elif intersection[0] == 2:
        #identical spheres
        #therefore no lower limit on the angle
        min_angle = 0.

    elif intersection[0] == 0:
        # Intersection is a point 
        min_angle = pi

    elif intersection[0] == 1:
        #print >>sys.stderr, "plane:", intersection

        o = array(intersection[1])
        n = array(intersection[2])
        r = intersection[3]

        p = get_point_on_circle(o, n, r)

        curr_pos = normalize(curr_pos)
        p = normalize(p)

        #print >>sys.stderr, "curr_pos:", curr_pos
        #print >>sys.stderr, "p:", p

        min_angle = acos(dot(curr_pos, p))


    return min_angle

def generate_loop(lengths):
    if len(lengths) < 3:
        print >> sys.stderr, "Not enough segments."
        sys.exit(1)

    i = 0

    prev_x = 0.
    prev_y = 0.
    prev_z = 0.

    origin = array([0., 0., 0,])
    print "lengths:", lengths

    while (i < len(lengths) - 2):
        trd = sum(lengths[i+1:])
        lrd = max(lengths[i+1:])

        curr_pos = array([prev_x, prev_y, prev_z])
        curr_distance = vec_distance(curr_pos, origin)
        cd = curr_distance
        
        if lrd - (trd - lrd) > 0:
            mind = lrd - (trd - lrd)
        else:
            mind = 0

        maxd = trd
        
        print >>sys.stderr, "i: %d cd: %f lrd: %f trd: %f mind: %f maxd: %f" % (i, cd, lrd, trd, mind, maxd)
        print >>sys.stderr, "lengths[i]: %d, cd: %f, mind: %f" % (lengths[i], cd, mind)

        if mind > abs(lengths[i] - cd):
            max_angle = pi - get_triangle_angle(lengths[i], cd, mind)
        else:
            max_angle = pi 

        if maxd >= abs(lengths[i] + cd):
            min_angle = 0
        else:
            min_angle = pi - get_triangle_angle(lengths[i], cd, maxd)

        print >>sys.stderr, "min_angle: %f max_angle: %f" % (min_angle, max_angle)
        #print >>sys.stderr, "lengths:", lengths

        #r = get_random_spherical_by_rotation(curr_pos, lengths[i], min_angle + 0.05, max_angle - 0.05)
        r = get_random_spherical_by_rotation(curr_pos, lengths[i], min_angle, max_angle)
        #r = get_random_spherical_by_rotation(curr_pos, lengths[i], 0., pi)
        #r = get_random_spherical(lengths[i])

        new_x = prev_x + r[0]
        new_y = prev_y + r[1]
        new_z = prev_z + r[2]

        print_segment((prev_x, prev_y, prev_z), (new_x, new_y, new_z))

        prev_x = new_x
        prev_y = new_y
        prev_z = new_z
            
        i += 1

    #connect_to_start(prev_x, prev_y, prev_z, lengths[-2:])
    curr_distance = vec_distance((prev_x, prev_y, prev_z), origin)
    print >>sys.stderr, "curr_distance:", curr_distance
    connect_in_3d((prev_x, prev_y, prev_z), origin, lengths[-2:])

def main():
    lengths = []

    #seed(1)

    for i in range(1, len(sys.argv)):
        lengths += [float(sys.argv[i])]

    segments = []
    #print_pymol_intro()
    generate_loop(lengths)

    #print_pymol_segments()
    #print_pymol_outro()
    #print_text()
    #print_angle_stats()

if __name__ == "__main__":
    main()
