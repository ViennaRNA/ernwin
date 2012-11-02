#!/usr/bin/python

import sys
import random as rand
import math
from optparse import OptionParser

import corgy.builder.models as cbm
import corgy.builder.reconstructor as cbr
import corgy.builder.stats as cbs
import corgy.graph.bulge_graph as cgb
import corgy.utilities.debug as cud
import corgy.utilities.vector as cuv

def get_random_stem_stats():
    '''
    Get a random stem from the statistics and return it.

    @return: A StemModel corresponding to a random stem
    '''
    # pick a length for the stem
    stem_stats = cbs.get_stem_stats()

    l = rand.choice(stem_stats.keys())
    return rand.choice(stem_stats[l])

def get_random_angle_stat():
    '''
    Create a random angle stats. This refers to the orienation
    of one helix with respect to another.

    @return: A random AngleStat
    '''
    return cbs.AngleStat('', 0, 0,
                             rand.uniform(0, math.pi),
                             rand.uniform(0, 2 * math.pi),
                             rand.uniform(0, 2 * math.pi),
                             rand.uniform(0, 100.),
                             rand.uniform(0, math.pi),
                             rand.uniform(0, 2 * math.pi))

def main():
    parser = OptionParser()

    parser.add_option('-n', '--nucleotides', dest='num_nucleotides', 
                      default=3, help="The number of nucleotides to sample", 
                      type='int')

    (options, args) = parser.parse_args()

    s1 = get_random_stem_stats()
    s2 = get_random_stem_stats()

    ang_stat = get_random_angle_stat()

    junction_len = 3
    s1_len = s1.define[1] - s1.define[0]
    s2_len = s2.define[1] - s2.define[0]

    # the new defines for the stem
    new_d1 = range(4)
    new_d2 = range(4)

    #Start with a little loop so that the traversal can get started
    loop_d = [0, 3]

    # shift the define of the first stem so that it is at the beginning
    new_d1[0] = 3
    new_d1[1] = 3 + s1_len

    # no residues between the two strands of stem
    new_d1[2] = s1_len + 1
    new_d1[3] = 2 * s1_len + 1

    # add the junction between the two stems
    new_d2[0] = new_d1[3] + junction_len + 1
    new_d2[1] = new_d2[0] + s2_len

    new_d2[2] = new_d2[1] + 1
    new_d2[3] = new_d2[2] + s2_len

    bg = cgb.BulgeGraph()

    bg.defines['s1'] = new_d1
    bg.defines['s2'] = new_d2
    bg.defines['b1'] = [new_d1[3], new_d2[0]]
    bg.defines['b2'] = loop_d

    bg.edges['s1'] = set(['b1', 'b2'])
    bg.edges['s2'] = set(['b1'])
    bg.edges['b1'] = set(['s1', 's2'])
    bg.edges['b2'] = set(['s1'])

    bg.weights['s1'] = 2
    bg.weights['s2'] = 2
    bg.weights['b1'] = 1
    bg.weights['b2'] = 1

    sm = cbm.SpatialModel(bg)

    sm.stem_defs['s1'] = s1
    sm.stem_defs['s2'] = s2
    sm.angle_defs['b1'] = ang_stat

    sm.traverse_and_build()

    sm.bg.output('out.bg')

    '''
    sm.add_stem('s1', s1, StemModel('start'), AngleStats(), (0,1))
    sm.add_stem('s2', s2, s1, ang_stats, (0,1))
    '''

if __name__ == "__main__":
    main()

