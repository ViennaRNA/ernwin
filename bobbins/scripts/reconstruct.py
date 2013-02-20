#!/usr/bin/python

import sys, pdb
import numpy as np
from optparse import OptionParser

import corgy.builder.reconstructor as rtor 
import corgy.builder.models as models

import corgy.utilities.debug as cud
import corgy.utilities.vector as cuv

from corgy.graph.bulge_graph import BulgeGraph

def main():
    parser = OptionParser()

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    parser.add_option('-l', '--loops', dest='loops', default=True, action='store_false', help='Toggle loop reconstruction')
    parser.add_option('','--drop-into-debugger', dest='drop_into_debugger', default=False, action='store_true')
    parser.add_option('-f', '--fragments', dest='fragments', default=False, action='store_true', help='Reconstruct using fragments.')
    parser.add_option('', '--output-file', dest='output_file', default='out.pdb', type='str')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        print >>sys.stderr, "Usage: ./reconstruct.py temp.comp"
        print >>sys.stderr, "Reconstruct a spatial model to full-atom accuracy."
        sys.exit(1)

    sm = models.SpatialModel(BulgeGraph(args[0]))
    sm.sample_native_stems()
    sm.create_native_stem_models()

    chain = rtor.reconstruct_stems(sm)

    if options.loops:
        if options.fragments:
            sm.sampled_from_bg()
            sampled_bulges = sm.get_sampled_bulges()

            for b in sm.bg.bulges():
                print >>sys.stderr, "reconstructing...", b
                if b not in sampled_bulges:
                    continue

                    print >> sys.stderr, "closed bulge:", b
                    (as1, as2) = sm.bg.get_bulge_angle_stats(b)
                    cud.pv('as1')
                    cud.pv('as2')

                    bulge_vec = np.array(cuv.spherical_polar_to_cartesian((as1.r1, as1.u1, as1.v1)))
                    #bulge_vec = np.array(cuv.spherical_polar_to_cartesian((as2.r1, as2.u1, as2.v1)))
                    size = sm.bg.get_bulge_dimensions(b)
                    cud.pv('(size, bulge_vec)')

                    best_bv_dist = 1000000.
                    best_bv = None
                    ang_type = sm.bg.get_angle_types(b)[0]

                    for ang_s in sm.angle_stats[size[0]][size[1]][ang_type]:
                        pot_bulge_vec = np.array(cuv.spherical_polar_to_cartesian((ang_s.r1, ang_s.u1, ang_s.v1)))
                        pot_bv_dist = cuv.magnitude(bulge_vec - pot_bulge_vec)
                        if pot_bv_dist < best_bv_dist:
                            best_bv_dist = pot_bv_dist
                            best_bv_ang_s = ang_s
                            best_bv = pot_bulge_vec

                    cud.pv('best_bv_dist')
                    cud.pv('best_bv')
                    cud.pv('bulge_vec')
                    cud.pv('best_bv_ang_s')
                    sm.angle_defs[b][ang_type] = best_bv_ang_s

                rtor.reconstruct_bulge_with_fragment(chain, sm, b)
            for l in sm.bg.loops():
                rtor.reconstruct_loop_with_fragment(chain, sm, l)
            for f in sm.bg.fiveprime():
                rtor.reconstruct_fiveprime_with_fragment(chain, sm, f)
            for t in sm.bg.threeprime():
                rtor.reconstruct_threeprime_with_fragment(chain, sm, t)
        else:
            try:
                #rtor.reconstruct_loop(chain, sm, 'b17', 0)
                rtor.reconstruct_loops(chain, sm, samples=40, consider_contacts=True)
            except Exception as e:
                if options.drop_into_debugger:
                    pdb.post_mortem()
                else:
                    raise

        #rtor.reconstruct_loops(chain, sm)

    rtor.replace_bases(chain, sm.bg.seq)
    rtor.output_chain(chain, options.output_file)

if __name__ == '__main__':
    main()

