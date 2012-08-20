import unittest

from numpy import pi, array, cos, sqrt
from numpy import cross, dot, allclose

from corgy.visual.pymol import PymolPrinter
from corgy.builder.config import Configuration

from corgy.utilities.vector import get_non_colinear_unit_vector, normalize
from corgy.utilities.vector import rotation_matrix, vec_angle, magnitude

from random import uniform

import os, sys


def third_side(b, c, A):
    '''
    Calculate the length of the third side of a triangle
    given two sides and the angle between them.

    For more information see:
    http://www.teacherschoice.com.au/maths_library/trigonometry/solve_trig_sas.htm
    
    @param b: The length of the first side.
    @param c: The length of the second side.
    @param A: The angle between the two sides.
    '''

    a = sqrt(b*b + c*c + 2. * b * c * cos(A))
    return a

def min_angles(lengths, angles):
    min_angle1 = angles[-1]
    max_angle1 = angles[-1]

    min_distance1 = third_side(lengths[-1], lengths[-2], angles[-1])
    max_distance1 = third_side(lengths[-1], lengths[-2], angles[-1])

    print "min_distance1:", min_distance1

def get_new_direction(prev_vec, length, angle):
    '''
    Get a new direction vector that is at a particular angle to the first one.

    @param prev_vec: The previous direction.
    @param length: The length of the vector.
    @param angle: The angle between the new vector and the previous one.
    '''
    vec2 = get_non_colinear_unit_vector(prev_vec)
    rot_vec = cross(prev_vec, vec2)

    rot1 = rotation_matrix(rot_vec, angle)
    rot2 = rotation_matrix(prev_vec, uniform(-pi, pi))

    rot_mat = dot(rot2, rot1)
    new_vec = length * normalize(dot(rot_mat, prev_vec))


    return new_vec

def dp(expr):
    return "print \"%s\", %s" % (expr, expr)

def is_reachable(p, c, lengths, angles):
    '''
    Is the vector from (0,0,0) to c reachable from the point p
    using a certain number of finite length segments, given angles. 

    The first length will be the length of the first segment (i.e.
    the one that will be touching c) and the first angle is the 
    angle separating the first segment from c.

    The next segment (lengths[1]) will touch lengths[0] with an angle
    of angles[1]. And so on.

    @param p: The destination point.
    @param c: The starting point and segment
    @param lengths: The lengths of the segments that may be used to connect
                    p to c
    @param angles: The angles between the segments.
    '''
    if len(lengths) == 1:
        d = dot(p-c, c)

        exec(dp("p"))
        exec(dp("c"))
        exec(dp("d"))

        exec dp("magnitude(c)")
        exec dp("lengths[0]")
        exec dp("cos(angles[0])")
        
        e = magnitude(c) * lengths[0] * cos(angles[0])
        exec dp("e")

        if allclose(dot(p-c, c), magnitude(c) * lengths[0] * cos(angles[0])):
            return True
        else:
            return False
    else:
        return False

class TestLoops(unittest.TestCase):
    def test_a(self):
        lengths = [1,1,1,1]
        angles = [pi/2, pi/2, pi/2, pi/2]

        min_angles(lengths, angles)

    def create_path(self, lengths, angles, n, pp, color):
        prev_loc = array([lengths[0], 0., 0.])
        prev_dir = prev_loc
        
        for i in range(1, n):
            prev_dir = get_new_direction(prev_dir, lengths[i], angles[i-1])
            prev_loc = prev_loc + prev_dir

       
        self.f.write("%d %f\n" % (n, magnitude(prev_loc)))

        pp.add_sphere(prev_loc, color=color, width=0.1)

    def test_is_reachable(self):
        prev_loc = array([2., 0., 0.])
        lengths = [2.]
        angles = [pi/2]

        prev_dir = prev_loc
        new_dir = get_new_direction(prev_dir, lengths[0], angles[0])

        self.assertEquals(is_reachable(prev_loc + new_dir, prev_loc, lengths, angles), True)
        self.assertEquals(is_reachable(prev_loc + new_dir + array([1., 0., 0.]), prev_loc, lengths, angles), False)


    def test_reachable(self):
        levels = 7
        n = 100

        all_colors = ['white', 'red', 'blue', 'green', 'yellow', 'purple']

        lengths = [2 for i in range(levels + 1)]
        angles = [pi/3 for i in range(levels + 1)]
        colors = [all_colors[i % len(all_colors)] for i in range(levels+1)]

        pp = PymolPrinter()
        pp.add_sphere(array([0., 0., 0.]), color='white', width=0.3)

        self.f = open(os.path.join(Configuration.test_output_dir, 'lengths.csv'), 'w')

        self.create_path(lengths, angles, 1, pp, 'red')

        #for i in range(2, levels+1):
        for i in range(6, 7):
            for k in range(n):
                self.create_path(lengths, angles, i, pp, colors[i])

        self.f.close()

        pp.dump_pymol_file(os.path.join(Configuration.test_output_dir, "angles"))

