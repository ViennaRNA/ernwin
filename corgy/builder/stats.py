#!/usr/bin/python

from random import uniform
from numpy import allclose
from corgy.utilities.data_structures import DefaultDict
from corgy.builder.config import Configuration

avg_stem_bp_length = 2.24
avg_twist_rotation_per_bp = 360 / 11.

loop_lengths = [ 
        (0., 0.),
        ( 7.0 , 9.0 ), 
        ( 7.459 , 9.33 ), 
        ( 7.774 , 8.945 ), 
        ( 8.102 , 8.985 ), 
        ( 6.771 , 8.182 ), 
        ( 6.465 , 7.533 ), 
        ( 6.435 , 7.676 ), 
        ( 6.605 , 8.987 ), 
        ( 8.396 , 9.367 ), 
        ( 12.13 , 18.68 ), 
        ( 19.76 , 22.32 ), 
        ( 11.57 , 14.59 ), 
        ( 8.702 , 8.744 ), 
        ( 15.46 , 15.46 ), 
        ( 15.0 , 30.0 ), 
        ( 15.0 , 30.0 ), 
        ( 15.0 , 30.0 ), 
        ( 15. , 30. ), 
        ( 15. , 30. ), 
        ( 15. , 30. ), 
        ( 15. , 30. ), 
        ( 15. , 30. ), 
        ( 15. , 30. ), 
        ( 33.02 , 33.02 ) ]

def get_loop_length(bg, key):
    if int(bg.defines[key][0]) > int(bg.defines[key][1]):
        loop_length = 1. 
    else:
        loop_length = int(bg.defines[key][1]) - int(bg.defines[key][0])

    return uniform(loop_lengths[loop_length][0], loop_lengths[loop_length][1])


class LoopStat:
    '''
    Class for storing the individual statistics about loops.

    phys_length: The length between the start and the centroid of the loop.
    '''
    def __init__(self, line=''):
        self.pdb_name = ''
        
        self.bp_length = 0
        self.phys_length = 0.

        if len(line) > 0:
            self.parse_line(line)

    def parse_line(self, line):
        '''
        Parse a line containing statistics about the shape of a stem.

        @param line: The line from the statistics file.
        '''

        parts =  line.strip().split()

        self.pdb_name = parts[1]

        self.bp_length = int(parts[2])
        self.phys_length = float(parts[3])

    def __str__(self):
        return "pdb_name: %s bp: %d phys_length: %f" % (self.pdb_name, self.bp_length, self.phys_length)

class StemStat:
    '''
    Class for storing the individual statistics about helices.

    Each stem will be defined by its base pair length. Two
    parameters are associated with each base pair length:

    phys_length: The physical length of such a helix
    twist_angle: The angle between its two twist segments
    '''
    def __init__(self, line=''):
        self.pdb_name = ''
        
        self.bp_length = 0
        self.phys_length = 0.

        self.twist_angle = 0.

        if len(line) > 0:
            self.parse_line(line)

    def parse_line(self, line):
        '''
        Parse a line containing statistics about the shape of a stem.

        @param line: The line from the statistics file.
        '''

        parts =  line.strip().split()

        self.pdb_name = parts[1]

        self.bp_length = int(parts[2])
        self.phys_length = float(parts[3])
        self.twist_angle = float(parts[4])

    def __str__(self):
        return "pdb_name: %s bp: %d phys_length: %f twist_angle: %f" % (self.pdb_name, self.bp_length, self.phys_length, self.twist_angle)

class AngleStat:
    '''
    Class for storing an individual statistic about inter-helical angles.
    '''

    def __init__(self, pdb_name='', dim1=0, dim2=0, u=0, v=0, t=0, r1=0, u1=0, v1=0):
        self.pdb_name = pdb_name
        self.dim1 = dim1
        self.dim2 = dim2

        self.u = u
        self.v = v
        self.t = t

        self.r1 = r1
        self.u1 = u1
        self.v1 = v1

    def __hash__(self):
        return id(self)

    def __eq__(self, a_s):
        '''
        Is this AngleStat equal to another one.
        '''
        if self.dim1 != a_s.dim1:
            return False
        if self.dim2 != a_s.dim2:
            return False
        if not allclose(self.u, a_s.u):
            return False
        if not allclose(self.v, a_s.v):
            return False
        if not allclose(self.t, a_s.t):
            return False
        if not allclose(self.r1, a_s.r1):
            return False
        if not allclose(self.u1, a_s.u1):
            return False
        if not allclose(self.v1, a_s.v1):
            return False

        return True

    def parse_line(self, line):
        parts = line.strip().split(' ')

        self.pdb_name = parts[1]

        self.dim1 = int(parts[2])
        self.dim2 = int(parts[3])

        self.u = float(parts[4])
        self.v = float(parts[5])
        self.t = float(parts[6])

        self.r1 = float(parts[7])
        self.u1 = float(parts[8])
        self.v1 = float(parts[9])

    def orientation_params(self):
        '''
        Return a tuple containing the parameters which specify the orientation
        of one stem with respect to another.

        @return: (u, v)
        '''

        return (self.u, self.v)

    def twist_params(self):
        '''
        Returns a tuple containing the parameters which specify the difference in
        twist between two stems.

        @return (u, v, t)
        '''
        return (self.u, self.v, self.t)

    def position_params(self):
        '''
        Return a tuple containing the parameters which specify the position
        of one stem with respect to another.

        @return (r1, u1, v1)
        '''
        return (self.r1, self.u1, self.v1)

    def __str__(self):
        str0 = "d1: %d d2: %d " % (self.dim1, self.dim2)
        str1 = "u: %f v: %f t: %f " % (self.u, self.v, self.t)
        str2 = "r1: %f u1: %f v1: %f" % (self.r1, self.u1, self.v1)
        return str0 + str1 + str2

class ConstructionStats:
    angle_stats = None
    stem_stats = None
    loop_stats = None

def get_angle_stats(filename=Configuration.stats_file):
    '''
    Load the statistics about inter the helix-helix orientations from a file.

    The file format should be as follows:

    angle pdb_name dim1 dim2 r u v t r1 u1 v1

    Where the parameters are as follows:

    angle: identifier for a type of statistics... should always just be 'angle'
    pdb_name: the name of the pdb file these statistics came from
    dim1: the smaller dimension of the bulge
    dim2: the larger dimension of the bulge
    r: the length of the second stem
    u: the polar angle of the orientation of the 2nd stem
    v: the azimuth of the orientation of the 2nd stem
    t: the orientation of the twist of the second stem
    
    r1: the distance of the start of the 2nd stem helix from the end of the 1st
        stem helix
    u1: the polar angle of the separation vector of the two helices
    v1: the azimuth of the separation vector of the two helices

    The azimuth is always defined with respect to the coordinate system defined
    by the stem1 helix axis vector and it's twist vector (the one adjacent to the
    bulge element).
    '''
    if ConstructionStats.angle_stats != None:
        return ConstructionStats.angle_stats

    ConstructionStats.angle_stats = DefaultDict(DefaultDict([]))

    f = open(filename, 'r')

    for line in f:
        if line.strip().find('angle') == 0:
            angle_stat = AngleStat()
            angle_stat.parse_line(line)
            ConstructionStats.angle_stats[angle_stat.dim1][angle_stat.dim2] += [angle_stat]

    f.close()

    return ConstructionStats.angle_stats


def get_stem_stats(filename=Configuration.stats_file):
    '''
    Load the statistics from the file.

    format:

    stem pdb_name bp_length phys_length twist_angle

    @param filename: The name of the file.
    '''

    if ConstructionStats.stem_stats != None:
        return ConstructionStats.stem_stats

    ConstructionStats.stem_stats = DefaultDict([])

    f = open(filename, 'r')

    for line in f:
        if line.strip().find('stem') == 0:
            stem_stat = StemStat(line)
            ConstructionStats.stem_stats[stem_stat.bp_length] += [stem_stat]

    f.close()

    return ConstructionStats.stem_stats


def get_loop_stats(filename=Configuration.stats_file):
    '''
    Load the statistics from the file.

    format:

    loop pdb_name bp_length phys_length

    @param filename: The name of the file.
    '''
    if ConstructionStats.loop_stats != None:
        return ConstructionStats.loop_stats

    ConstructionStats.loop_stats = DefaultDict([])

    f = open(filename, 'r')

    for line in f:
        if line.strip().find('loop') == 0:
            loop_stat = LoopStat(line)
            ConstructionStats.loop_stats[loop_stat.bp_length] += [loop_stat]

    f.close()

    return ConstructionStats.loop_stats
