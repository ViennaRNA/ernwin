import unittest

from numpy import pi, array, cos, sqrt, linspace
from numpy import cross, dot, allclose

import numpy as np

from math import acos
from copy import copy

from borgy.visual.pymol import PymolPrinter
from borgy.builder.config import Configuration

from borgy.utilities.vector import get_non_colinear_unit_vector, normalize
from borgy.utilities.vector import rotation_matrix, vec_angle, magnitude
from borgy.utilities.vector import get_random_vector

from random import uniform, random

import os, sys

import matplotlib.pyplot as plt


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

def cos_angle_from_sides(aside1, aside2, oside):
    return (aside1 ** 2 + aside2 ** 2 - oside ** 2) / (2 * aside1 * aside2)

def oside_from_sides_and_angle(aside1, aside2, angle):
    return sqrt(aside1**2 + aside2**2 - 2. * aside1 * aside2 * cos(angle))

def min_angles(lengths, angles):
    min_angle1 = angles[-1]
    max_angle1 = angles[-1]

    min_distance1 = third_side(lengths[-1], lengths[-2], angles[-1])
    max_distance1 = third_side(lengths[-1], lengths[-2], angles[-1])

def get_new_direction(prev_vec, length, angle, rot_angle=None, pp_vec=None):
    '''
    Get a new direction vector that is at a particular angle to the first one.

    @param prev_vec: The previous direction.
    @param length: The length of the vector.
    @param angle: The angle between the new vector and the previous one.
    '''
    if rot_angle == None:
        rot_angle = uniform(-pi, pi)

    if pp_vec == None:
        vec2 = get_non_colinear_unit_vector(prev_vec)
    else:
        vec2 = pp_vec


    rot_vec = cross(prev_vec, vec2)

    rot1 = rotation_matrix(rot_vec, angle)
    rot2 = rotation_matrix(prev_vec, rot_angle)

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
        
        e = magnitude(c) * lengths[0] * cos(angles[0])

        if allclose(dot(p-c, c), magnitude(c) * lengths[0] * cos(angles[0])):
            return True
        else:
            return False

    elif len(lengths) == 2:
        len_p1_c = oside_from_sides_and_angle(lengths[0], lengths[1], pi - angles[1])
        cos_ax = cos_angle_from_sides(len_p1_c, lengths[0], lengths[1])

        a_pc = vec_angle(p, c)
        
        if cos_ax > 1. or cos_ax < -1.:
            return False

        a_x = acos(cos_ax)

        if allclose(magnitude(p - c), len_p1_c):
            return True
        return False

    #print "len(lengths):", len(lengths)
    if len(lengths) >= 3:
        len_p1_c = oside_from_sides_and_angle(lengths[1], lengths[2], pi - angles[2])
        cos_ax = cos_angle_from_sides(len_p1_c, lengths[1], lengths[2])

        if cos_ax > 1. or cos_ax < -1.:
            return False

        a_x = acos(cos_ax)

        len_p1_min = oside_from_sides_and_angle(lengths[0], len_p1_c, pi - angles[1] - a_x)
        len_p1_max = oside_from_sides_and_angle(lengths[0], len_p1_c, pi - angles[1] + a_x)

        #len_p1_min = sqrt(lengths[0] ** 2 + len_p1_c ** 2 - 2 * lengths[0] * len_p1_c * cos(pi - angles[1] - a_x))
        #len_p1_max = sqrt(lengths[0] ** 2 + len_p1_c ** 2 - 2 * lengths[0] * len_p1_c * cos(pi - angles[1] + a_x))

        cos_a_p1_max = cos_angle_from_sides(len_p1_min, lengths[0], len_p1_c)
        cos_a_p1_min = cos_angle_from_sides(len_p1_max, lengths[0], len_p1_c)

        a_p1_max = acos(cos_a_p1_max)
        a_p1_min = acos(cos_a_p1_min)

        '''
        print "len_p1_min:", len_p1_min
        print "len_p1_max:", len_p1_max
        print "a_p1_min:", a_p1_min
        print "a_p1_max:", a_p1_max
        print "p:", magnitude(p)
        '''

        tol = 1e-7
        if len(lengths) == 3:
            if len_p1_min - tol <= magnitude(p) and magnitude(p) <= len_p1_max + tol:
                return True

            return False

    if len(lengths) >= 4:
        cos_ay = cos_angle_from_sides(lengths[2], len_p1_c, lengths[1])
        cos_az_min = cos_angle_from_sides(len_p1_c, len_p1_min, lengths[0])
        cos_az_max = cos_angle_from_sides(len_p1_c, len_p1_max, lengths[0])

        ay = acos(cos_ay)
        az_min = acos(cos_az_min)
        az_max = acos(cos_az_max)

        aw_min = pi - angles[3] - ay - az_min
        aw_max = pi - angles[3] + ay - az_max
        possible_aw_maxs = [pi - angles[3] + ay - az_max, pi - angles[3] - ay - az_max, pi-angles[3] + ay + az_max, pi - angles[3] - ay + az_max]

        len_p2_maxs = [oside_from_sides_and_angle(len_p1_max, lengths[3], aw) for aw in possible_aw_maxs]
        len_p2_max = max([oside_from_sides_and_angle(len_p1_max, lengths[3], aw) for aw in possible_aw_maxs])
            
        len_p2_min = oside_from_sides_and_angle(len_p1_min, lengths[3], aw_min)
        #len_p2_max = oside_from_sides_and_angle(len_p1_max, lengths[3], aw_max)

        #print "angles[3]:", pi - angles[3]
        print "az_min:", az_min, "az_max:", az_max, "ay:", ay
        #print "ay:", ay
        #print "max_angle:", acos(cos_angle_from_sides(len_p1_max, lengths[3], 3.371348))
        #print "other max_angle:", pi-angles[3]+ay+az_max
        #print "angles[3] + ay:", pi-angles[3] + az_max
        #print "----------------"
        #print "p1_min:", len_p1_min, "p1_max:", len_p1_max
        #print "p2_min:", len_p2_min, "p2_max:", len_p2_max
        #print "----------------"
        #print "magnitude:", magnitude(p), "lengths:", lengths, "angles:", [180. * a / pi for a in angles]

        if len(lengths) == 4:
            if len_p2_min <= magnitude(p) and magnitude(p) <= len_p2_max:
                return True

            return False

        pass


class TestLoops(unittest.TestCase):
    def setUp(self):
        np.seterr(all='ignore')

    def test_a(self):
        lengths = [1,1,1,1]
        angles = [pi/2, pi/2, pi/2, pi/2]

        min_angles(lengths, angles)

    def create_path(self, lengths, angles, n, pp, color, limits=[]):
        prev_loc = array([lengths[0], 0., 0.])
        prev_dir = prev_loc
        pp_dir = None
        
        for i in range(1, n):
            if i in limits:
                if random() < 0.5:
                    prev_dir = get_new_direction(prev_dir, lengths[i], angles[i-1], rot_angle=0, pp_vec = array([10., 11., 0.]))
                else:
                    prev_dir = get_new_direction(prev_dir, lengths[i], angles[i-1], rot_angle=0, pp_vec = array([10., -11., 7.]))

            else:
                prev_dir = get_new_direction(prev_dir, lengths[i], angles[i-1])

            prev_loc = prev_loc + prev_dir

       
        self.f.write("%d %f\n" % (n, magnitude(prev_loc)))

        pp.add_sphere(prev_loc, color=color, width=0.1)

    def plot_directions(self, ax, dir2, dir1, loc1, length, angle, label):
        i_s = []
        loc_s = []

        for i in linspace(-pi, pi, 50):
            dir3 = get_new_direction(dir2, length, angle, rot_angle=i, pp_vec=dir1)
            loc = loc1 + dir1 + dir2 + dir3
            i_s += [i]
            #i_s += [vec_angle(dir3, array([1, 0., 0.]))]
            loc_s += [magnitude(loc)]

        ax.plot(i_s, loc_s, label=label)

    def test_min_distance(self):
        lengths = [uniform(1,1) for i in xrange(10)]
        angles = [uniform(0, pi / 2.) for i in xrange(10)]
        angles = [82.5093958289398, 59.41108471194514, 37.123527549300654, 70.38114153481955]
        angles = [a * pi / 180. for a in angles]

        '''
        lengths = [uniform(1,1) for i in xrange(10)]
        angles = [pi/2. for i in xrange(10)]
        '''

        ax = plt.subplot(1,1,1)

        c = array([0., 0., 0.])
        dir0 = array([lengths[0], 0., 0.])

        dir1 = get_new_direction(dir0, lengths[1], angles[1])
        dir2 = get_new_direction(dir1, lengths[2], angles[2])
        dir3 = get_new_direction(dir2, lengths[3], angles[3])

        self.assertTrue(allclose(vec_angle(dir0, dir1), angles[1]))
        self.assertTrue(allclose(vec_angle(dir1, dir2), angles[2]))
        self.assertTrue(allclose(vec_angle(dir2, dir3), angles[3]))

        self.plot_directions(ax, dir1, dir0, c, lengths[2], angles[2], "dir2")
        self.plot_directions(ax, dir2, dir1, dir0, lengths[3], angles[3], "dir3")
        #self.plot_directions(ax, dir3, dir2, dir1, lengths[3], angles[3], "dir3")

        self.assertEquals(is_reachable(c + dir0 + dir1, c, lengths[:2], angles[:2]), True)
        self.assertEquals(is_reachable(c + dir0 + dir1 + dir2, c + dir0, lengths[1:3], angles[1:3]), True)

        min_lens = [1000. for i in range(10)]
        max_lens = [0. for i in range(10)]

        min_angles = [1000. for i in range(10)]
        max_angles = [0. for i in range(10)]

        self.assertEquals(is_reachable(c + dir0 + dir1 + dir2, c, lengths[:3], angles[:3]), True)
        self.assertEquals( is_reachable(c + dir0 + dir1 + dir2 + dir3, c, lengths[:4], angles[:4]), True)

        for k in range(200):
            dirs = [dir0]
            lens = [magnitude(c + dir0)]
            p_angles = [0.]

            for i in range(1, 4):
                angle1 = uniform(-pi, pi)
                angle2 = uniform(-pi, pi)

                if i == 3:
                    dirs += [get_new_direction(dirs[i-1], lengths[i], angles[i], rot_angle=angle1, pp_vec = dirs[i-2])]
                elif i == 2:
                    dirs += [get_new_direction(dirs[i-1], lengths[i], angles[i], rot_angle=angle2, pp_vec = dirs[i-2])]
                else:
                    dirs += [get_new_direction(dirs[i-1], lengths[i], angles[i])]


                loc = copy(c)
                
                for j in range(0, i+1):
                    loc += dirs[j]

                #self.assertEquals(is_reachable(loc, c, lengths[:i+1], angles[:i+1]), True)

                #print "c:", c
                lens += [magnitude(loc)]
                p_angles += [vec_angle(dir0, loc)]

                if lens[i] > max_lens[i]:
                    if i == 3:
                        print "angle1:", angle1, "angle2:", angle2
                        is_reachable(loc, c, lengths[:i+1], angles[:i+1])
                    max_lens[i] = lens[i]
                if lens[i] < min_lens[i]:
                    min_lens[i] = lens[i]

                if p_angles[i] > max_angles[i]:
                    max_angles[i] = p_angles[i]
                if p_angles[i] < min_angles[i]:
                    min_angles[i] = p_angles[i]

        for i in range(1, 4):
            print "min_len %d: %f max_len: %f" % (i-1, min_lens[i], max_lens[i])
            #print "min_angle %d: %f max_angle: %f" % (i, min_angles[i], max_angles[i])

        #self.plot_directions(ax, dir3, dir2, dir0 + dir1, lengths[3], angles[3], "dir3")
        print "angles:", angles

        handles, labels = ax.get_legend_handles_labels()
        ax.legend(handles, labels)

        #plt.show()
        

    def test_is_reachable(self):
        prev_loc = array([2., 0., 0.])
        lengths = [2., 2., 2.]
        angles = [pi/2, pi/2, pi/2]

        prev_dir = prev_loc
        new_dir = get_new_direction(prev_dir, lengths[0], angles[0])

        self.assertEquals(is_reachable(prev_loc + new_dir, prev_loc, lengths[:1], angles[:1]), True)
        self.assertEquals(is_reachable(prev_loc + new_dir + array([1., 0., 0.]), prev_loc, lengths[:1], angles[:1]), False)

        new_dir1 = get_new_direction(new_dir, lengths[1], angles[1])
        self.assertEquals(is_reachable(prev_loc + new_dir + new_dir1, prev_loc, lengths[:2], angles[:2]), True)
        self.assertEquals(is_reachable(prev_loc + new_dir + get_random_vector(), prev_loc, lengths[:2], angles[:2]), False)


    def test_reachable(self):
        levels = 7
        n = 400

        all_colors = ['white', 'red', 'blue', 'green', 'yellow', 'purple']

        lengths = [2 for i in range(levels + 1)]
        angles = [pi/2. for i in range(levels + 1)]

        #lengths = [(i+1)*2 for i in range(levels + 1)]
        #angles = [(i * uniform(.8, 1.) * pi/3) for i in range(levels + 1)]

        colors = [all_colors[i % len(all_colors)] for i in range(levels+1)]

        pp = PymolPrinter()
        pp.add_sphere(array([0., 0., 0.]), color='white', width=0.3)

        self.f = open(os.path.join(Configuration.test_output_dir, 'lengths.csv'), 'w')

        self.create_path(lengths, angles, 1, pp, 'red')

        #for i in range(2, levels+1):
        xmin=4
        xmax=5

        for i in range(xmin, xmax):
            for k in range(n):
                self.create_path(lengths, angles, i, pp, colors[i], limits=[j for j in range(i-1)])

        for i in range(xmin-1, xmax-1):
            for k in range(n):
                self.create_path(lengths, angles, i, pp, colors[i])


        self.f.close()

        pp.dump_pymol_file(os.path.join(Configuration.test_output_dir, "angles"))

