#!/usr/bin/python

import sys, pdb
import numpy as np
import itertools as it
from optparse import OptionParser

import forgi.threedee.model.stats as cbs
import fess.builder.reconstructor as rtor
import fess.builder.reconstructor as rtor 
import fess.builder.models as models

import forgi.utilities.debug as fud
import forgi.threedee.utilities.vector as cuv

import forgi.threedee.model.coarse_grain as ftmc

def main():
    parser = OptionParser()

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    parser.add_option('-l', '--loops', dest='loops', default=True, action='store_false', help='Toggle loop reconstruction')
    parser.add_option('','--drop-into-debugger', dest='drop_into_debugger', default=False, action='store_true')
    parser.add_option('-f', '--fragments', dest='fragments', default=False, action='store_true', help='Reconstruct using fragments.')
    parser.add_option('', '--output-file', dest='output_file', default='out.pdb', type='str')
    parser.add_option('-s', '--samples', dest='samples', default=10, type='int', help='The number of samples to get from Barnacle')
    parser.add_option('-a', '--average_atoms', dest='average_atoms', default=False, action='store_true', help='Reconstruct using the positions of the average atoms')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        print >>sys.stderr, "Usage: ./reconstruct.py temp.comp"
        print >>sys.stderr, "Reconstruct a spatial model to full-atom accuracy."
        sys.exit(1)

    print >>sys.stderr, "reconstructing model:", args[0]
    sm = models.SpatialModel(ftmc.CoarseGrainRNA(args[0]))
    sm.sample_native_stems()
    sm.create_native_stem_models()

    if not options.average_atoms:
        chain = rtor.reconstruct_stems(sm)

    if options.loops:
        if options.average_atoms:
            chain = rtor.reconstruct_from_average(sm)
            rtor.output_chain(chain, 'temp.pdb')
            return

        elif options.fragments:
            sm.sampled_from_bg()
            sampled_bulges = sm.get_sampled_bulges()

            '''
            rtor.reconstruct_bulge_with_fragment(chain, sm, 'x2')
            rtor.output_chain(chain, options.output_file)
            sys.exit(1)
            '''

            for b in it.chain(sm.bg.iloop_iterator(), sm.bg.mloop_iterator()):

                if b not in sm.bg.sampled.keys():
                    print >>sys.stderr, "interpolating... %s" % (b)
                    (as1, as2) = sm.bg.get_bulge_angle_stats(b)

                    bulge_vec = np.array(cuv.spherical_polar_to_cartesian((as1.r1, as1.u1, as1.v1)))
                    #bulge_vec = np.array(cuv.spherical_polar_to_cartesian((as2.r1, as2.u1, as2.v1)))
                    size = sm.bg.get_bulge_dimensions(b)

                    best_bv_dist = 1000000.
                    best_bv = None
                    connections = list(sm.bg.edges[b])
                    ang_type = sm.bg.connection_type(b, connections)

                    for ang_s in cbs.get_angle_stats()[(size[0], size[1], ang_type)]:
                        pot_bulge_vec = np.array(cuv.spherical_polar_to_cartesian((ang_s.r1, ang_s.u1, ang_s.v1)))
                        pot_bv_dist = cuv.magnitude(bulge_vec - pot_bulge_vec)
                        if pot_bv_dist < best_bv_dist:
                            best_bv_dist = pot_bv_dist
                            best_bv_ang_s = ang_s
                            best_bv = pot_bulge_vec

                    sm.angle_defs[b][ang_type] = best_bv_ang_s
                    rtor.reconstruct_bulge_with_fragment(chain, sm, b, move_all_angles=True)
                    continue

                rtor.reconstruct_bulge_with_fragment(chain, sm, b)
            for l in sm.bg.hloop_iterator():
                rtor.reconstruct_loop_with_fragment(chain, sm, l)
            for f in sm.bg.floop_iterator():
                rtor.reconstruct_fiveprime_with_fragment(chain, sm, f)
            for t in sm.bg.tloop_iterator():
                rtor.reconstruct_threeprime_with_fragment(chain, sm, t)
        else:
            try:
                rtor.output_chain(chain, 'stems.pdb')
                #rtor.reconstruct_loop(chain, sm, 'b1', 0, samples=options.samples, consider_contacts=False)
                rtor.reconstruct_loops(chain, sm, samples=options.samples, consider_contacts=False)
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

