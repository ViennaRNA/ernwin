#!/usr/bin/python

from __future__ import absolute_import
from __future__ import print_function
import sys
from optparse import OptionParser

import itertools as it
import tess.threedee.model.coarse_grain as ttmc
import borgy.utilities.vector as cuv
import borgy.graph.graph_pdb as cgg
from six.moves import range

def main():
    usage = './virtual_res_distances temp.comp'
    usage += """
    List the distances between each virtual residue in the structure.
    """

    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    bg = ttmc.CoarseGrainRNA(args[0])
    bg.calc_bp_distances()
    stems = bg.stems()
    adjacent = 0

    mult = 7.

    for (s1, s2) in it.combinations(bg.stems(), 2):
        s1_len = bg.stem_length(s1)
        s2_len = bg.stem_length(s2)

        adjacent = 0

        if bg.are_any_adjacent_stems(s1, s2):
            adjacent = 1

        for k in range(s1_len):
            for l in range(s2_len):
                (v1_p, v1_v, v1_l, v1_r) = cgg.virtual_res_3d_pos(bg, s1, k)
                (v2_p, v2_v, v2_l, v2_r) = cgg.virtual_res_3d_pos(bg, s2, l)

                closest_points = cuv.line_segment_distance(v1_p, 
                                                           v1_p + mult * v1_v,
                                                           v2_p,
                                                           v2_p + mult * v2_v)

                closest_distance = cuv.magnitude(closest_points[1] - closest_points[0])
                vres_distance = bg.calc_vres_distance(s1, k, s2, l)

                print("dist:", cuv.magnitude((v2_p + mult * v2_v) - (v1_p + mult * v1_v)), closest_distance, adjacent, s1, k, s2, l, vres_distance)

            for l in range(k+1, s1_len):
                (v1_p, v1_v, v1_l, v1_r) = cgg.virtual_res_3d_pos(bg, s1, k)
                (v2_p, v2_v, v2_l, v2_r) = cgg.virtual_res_3d_pos(bg, s2, l)

                closest_points = cuv.line_segment_distance(v1_p, 
                                                           v1_p + mult * v1_v,
                                                           v2_p,
                                                           v2_p + mult * v2_v)

                closest_distance = cuv.magnitude(closest_points[1] - closest_points[0])
                vres_distance = bg.calc_vres_distance(s1, k, s2, l)

                print("internal_dist:", cuv.magnitude((v2_p + mult * v2_v) - (v1_p + mult * v1_v)), closest_distance, 1, s1, k, s1, l, vres_distance)

                        
if __name__ == '__main__':
    main()

