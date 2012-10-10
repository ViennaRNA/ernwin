#!/usr/bin/python

import sys, os
from optparse import OptionParser

import numpy as np

import corgy.graph.bulge_graph as cgb
import corgy.graph.graph_pdb as cgg
import corgy.builder.config as cbc
import corgy.visual.pymol as cvp
import corgy.utilities.vector as cuv

def main():
    usage = """
Visualize the data about the location of surrounding virutal
residues.

usage: %prog [options] data_file
"""
    parser = OptionParser(usage = usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        print >>sys.stderr, "Lack of data_file argument."
        parser.print_help()
        sys.exit(1)

    bg = cgb.BulgeGraph(os.path.join(cbc.Configuration.test_input_dir, "1y26/graph", "temp.comp"))

    stem = 's3' 
    stem_len = bg.defines[stem][1] - bg.defines[stem][0] + 1

    for i in xrange(stem_len):
        res_type = cgg.get_residue_type(i, stem_len)
        #print "res_type:", res_type

    pp = cvp.PymolPrinter()
    pp.add_stem(bg, stem)
    #pp.coordinates_to_pymol(bg)


    f = open(args[0], 'r')
    for line in f.readlines():
        parts = line.split(' ')
        res_type = int(parts[0])
        x = float(parts[2])
        y = float(parts[3])
        z = float(parts[4])

        for i in range(stem_len):
            s1_start = cgg.pos_to_spos(bg, stem, i, stem, 0)
            s1_end = cgg.pos_to_spos(bg, stem, i, stem, stem_len - 1)

            i_type = cgg.get_residue_type(i, stem_len)

            if i_type == res_type:
                pos = cgg.spos_to_pos(bg, stem, i, np.array([x,y,z]))
                #print "pos:", pos

                if x >= s1_start[0] and x <= s1_end[0]:
                    #pp.add_sphere(pos, "red", 1.7, '')
                    pp.add_sphere(pos, "red", 0.1, '')

    pp.output_pymol_file()

if __name__ == '__main__':
    main()

