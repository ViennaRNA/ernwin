#!/usr/bin/python

import sys
import itertools as it

from optparse import OptionParser

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.vector as cuv
import forgi.utilities.debug as fud

from collections import defaultdict

def output_long_range_distances(bg):
    for key1 in bg.longrange.keys():
        for key2 in bg.longrange[key1]:
            if bg.has_connection(key1, key2):
                continue

            #point1 = bg.get_point(key1)
            #point2 = bg.get_point(key2)

            (i1,i2) = cuv.line_segment_distance(bg.coords[key1][0], bg.coords[key1][1],
                                             bg.coords[key2][0], bg.coords[key2][1])

            vec1 = bg.coords[key1][1] - bg.coords[key2][0]
            basis = cuv.create_orthonormal_basis(vec1)
            coords2 = cuv.change_basis(i2 - i1, basis, cuv.standard_basis)
            (r, u, v) = cuv.spherical_polar_to_cartesian(coords2)

            
            seq1 = 'x'
            seq2 = 'x'

            '''
            if key1[0] != 's' and key[0] != 'i':
                seq1 = bg.get_seq(key1)
            if key2[0] != 's' and key2[0] != 'i':
                seq2 = bg.get_seq(key2)
            '''

            print "%s %s %d %s %s %d %f %s %s %s %f" % (key1, 
                                         key1[0], 
                                         bg.get_length(key1),
                                         key2, 
                                         key2[0],
                                         bg.get_length(key2),
                                         cuv.magnitude(i2-i1),
                                         seq1, seq2, "Y", v)

def output_all_distances(bg):
    for (key1, key2) in it.permutations(bg.defines.keys(), 2):
        if bg.has_connection(key1, key2):
            continue

        longrange = "N"

        if key2 in bg.longrange[key1]:
            longrange = "Y"

        #point1 = bg.get_point(key1)
        #point2 = bg.get_point(key2)

        try:
            (i1,i2) = cuv.line_segment_distance(bg.coords[key1][0], bg.coords[key1][1],
                                             bg.coords[key2][0], bg.coords[key2][1])


            if abs(cuv.magnitude(i2 - i1)) < 0.000001:
                continue

            vec1 = bg.coords[key1][1] - bg.coords[key1][0]
            '''
            basis = cuv.create_orthonormal_basis(vec1)
            coords2 = cuv.change_basis(i2 - i1, basis, cuv.standard_basis)
            (r, u, v) = cuv.spherical_cartesian_to_polar(coords2)
            '''
            v = cuv.vec_angle(vec1, i2 - i1)

        except KeyError as ke:
            #print >>sys.stderr, 'Skipping %s or %s.' % (key1, key2)
            continue

        seq1 = 'x'
        seq2 = 'x'


        '''
        receptor_angle = 0.
        if bg.get_type(key1) != 's' and bg.get_type(key1) != 'i' and bg.get_length(key1) > 1:
            seq1 = bg.get_seq(key1)
        if bg.get_type(key2) != 's' and bg.get_type(key2) != 'i'and bg.get_length(key2) > 1:
            seq2 = bg.get_seq(key2)
        if bg.get_type(key1) == 'l' and bg.get_type(key2) == 's':
            receptor_angle = cgg.receptor_angle(bg, key1, key2)
        '''

        print "%s %s %d %s %s %d %f %s %s %s %f" % (key1, 
                                     key1[0], 
                                     bg.get_length(key1),
                                     key2, 
                                     key2[0],
                                     bg.get_length(key2),
                                     cuv.magnitude(i2-i1),
                                     seq1, seq2, longrange, v)
def main():
    parser = OptionParser()

    parser.add_option('-a', '--all-pairs', dest='all_pairs', default=False, action='store_true', help='Print all interactions')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        print >>sys.stderr, "Usage: ./long_range_distances.py temp.comp"
        sys.exit(1)

    bg = ftmc.CoarseGrainRNA(args[0])

    if options.all_pairs:
        output_all_distances(bg)
    else: 
        output_long_range_distances(bg)

    

if __name__ == '__main__':
    main()

