#!/usr/bin/python

import sys
import random as rand
import math
from optparse import OptionParser

import corgy.builder.stats as cbs
import corgy.utilities.debug as cud

def get_random_stem():
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

    parser.add_option('-n', '--nucleotides', dest='num_nucleotides', default=3, help="The number of nucleotides to sample", type='int')

    (options, args) = parser.parse_args()

    s1 = get_random_stem()
    s2 = get_random_stem()

    ang_stat = get_random_angle_stat()

if __name__ == "__main__":
    main()

