#!/usr/bin/python

import sys
from optparse import OptionParser

import corgy.builder.reconstructor as rtor 
import corgy.builder.models as models

from corgy.graph.bulge_graph import BulgeGraph

def main():
    parser = OptionParser()

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    parser.add_option('-l', '--loops', dest='loops', default=True, action='store_false', help='Toggle loop reconstruction')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        print >>sys.stderr, "Usage: ./reconstruct.py temp.comp"
        print >>sys.stderr, "Reconstruct a spatial model to full-atom accuracy."
        sys.exit(1)

    sm = models.SpatialModel(BulgeGraph(args[0]))
    sm.sample_native_stems()
    sm.create_native_stem_models()

    chain = rtor.reconstruct_stems(sm)
    rtor.replace_bases(chain, sm.bg.seq)

    if options.loops:
        rtor.reconstruct_loops(chain, sm)

    rtor.output_chain(chain, 'out.pdb')

if __name__ == '__main__':
    main()

