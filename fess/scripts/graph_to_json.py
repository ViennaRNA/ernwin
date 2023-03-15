#!/usr/bin/python

from __future__ import absolute_import
from __future__ import print_function
import sys

import tess.threedee.model.coarse_grain as ttmc
import borgy.visual.force_graph as cvf

from optparse import OptionParser

def main():
    parser = OptionParser()

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        print("./graph_to_json temp.graph", file=sys.stderr)
        sys.exit(1)

    bg = ttmc.CoarseGrainRNA(args[0])
    json_str = cvf.convert_graph_to_fancy_json(bg)
    print(json_str)

if __name__ == '__main__':
    main()

