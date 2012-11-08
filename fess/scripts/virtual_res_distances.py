#!/usr/bin/python

import sys
from optparse import OptionParser

import corgy.graph.bulge_graph as cgb
import corgy.utilities.vector as cuv
import corgy.graph.graph_pdb as cgg

def main():
    usage = './virtual_res_distances temp.comp'
    usage = """
    List the distances between each virtual residue in the structure.
    """

    parser = OptionParser()

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    bg = cgb.BulgeGraph(args[0])
    #stems = [d for d in bg.defines.keys() if d[0] == 's']
    stems = [d for d in bg.defines.keys() if (bg.weights[d] == 2 or bg.weights[d] == 0)]

    mult = 7.

    for i in range(len(stems)):
        #s1_len = bg.defines[stems[i]][1] - bg.defines[stems[i]][0] + 1
        s1_len = bg.stem_length(stems[i])
        for j in range(i + 1, len(stems)):
            #s2_len = bg.defines[stems[j]][1] - bg.defines[stems[j]][0] + 1
            s2_len = bg.stem_length(stems[i])
            for k in range(s1_len):
                for l in range(s2_len):
                    (v1_p, v1_v) = cgg.virtual_res_3d_pos(bg, stems[i], k)
                    (v2_p, v2_v) = cgg.virtual_res_3d_pos(bg, stems[j], l)

                    print "dist:", cuv.magnitude((v2_p + mult * v2_v) - (v1_p + mult * v1_v)), stems[i], k, stems[j], l

                for l in range(k+1, s1_len):
                    (v1_p, v1_v) = cgg.virtual_res_3d_pos(bg, stems[i], k)
                    (v2_p, v2_v) = cgg.virtual_res_3d_pos(bg, stems[i], l)

                    print "internal_dist:", cuv.magnitude((v2_p + mult * v2_v) - (v1_p + mult * v1_v)), stems[i], k, stems[i], l

                        
if __name__ == '__main__':
    main()

