#!/usr/bin/python

import sys, os
from optparse import OptionParser

import numpy as np

import corgy.graph.bulge_graph as cgb
import corgy.graph.graph_pdb as cgg

import corgy.builder.config as cbc

import corgy.utilities.vector as cuv

def connected_stems(bg, s1, s2):
    for edge in bg.edges[s1]:
        if s2 in bg.edges[edge]:
            return True
    return False

def main():
    usage = """
Collect the information about how stems are surrounded by nucleotides.

For each nucleotide in each stem, the position 
(relative to that nucleotide) of each surrounding
nucleotide (which is also part of a stem) will be recorded. 

The information will be classified according to the type
of residue that the base is.

usage: %prog [options] temp.comp
"""
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        print >>sys.stderr, "Missing graph file."
        parser.print_help()
        sys.exit(1)

    bg = cgb.BulgeGraph(args[0])

    #stems = [d for d in bg.defines.keys() if d[0] == 's']
    stems = [d for d in bg.defines.keys() if (bg.weights[d] == 2 or bg.weights[d] == 0)]

    for i in range(len(stems)):
        if stems[i][0] != 's':
            # the end nucleotides of bulges overlap those of stems
            # so their length needs to be one less than that of stems
            s1_len = bg.defines[stems[i]][1] - bg.defines[stems[i]][0]
            k_start = 1
        else:
            s1_len = bg.defines[stems[i]][1] - bg.defines[stems[i]][0] + 1
            k_start = 0

        for j in range(len(stems)):
            if i == j:
                continue

            #if connected_stems(bg, stems[i], stems[j]):
            #    continue

            if stems[j][0] != 's':
                # see comment for s1_len above
                s2_len = bg.defines[stems[j]][1] - bg.defines[stems[j]][0]
                l_start = 1
            else:
                s2_len = bg.defines[stems[j]][1] - bg.defines[stems[j]][0] + 1
                l_start = 0

            for k in range(k_start, s1_len):
                s1_start = cgg.pos_to_spos(bg, stems[i], k, stems[i], 0)
                s1_end = cgg.pos_to_spos(bg, stems[i], k, stems[i], s1_len - 1)

                '''
                print >>sys.stderr, stems[i], k, s1_start, s1_end
                print >>sys.stderr, s1_len
                '''

                #print "s1_start:", s1_start, "s1_end:", s1_end
                for l in range(l_start, s2_len):
                    r1_type = cgg.get_residue_type(k, s1_len)
                    r2_spos = cgg.pos_to_spos(bg, stems[i], k, stems[j], l)

                    if cuv.magnitude(r2_spos) < 400. and r2_spos[0] > s1_start[0] and r2_spos[0] < s1_end[0]:
                        print r1_type, cuv.magnitude(r2_spos), " ".join(map(str, r2_spos)), bg.name, stems[i], k, stems[j], l


if __name__ == '__main__':
    main()

