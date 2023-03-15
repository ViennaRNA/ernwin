#!/usr/bin/python

from __future__ import absolute_import
from __future__ import print_function
import sys
from optparse import OptionParser

import fess.builder.energy as fbe
import forgi.threedee.model.coarse_grain as ftmc
import forgi.utilities.debug as cud
from six.moves import map
from six.moves import range

def main():
    usage = """
    python create_colinearity_measures.py temp.cg
    """
    num_args= 1
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')
    parser.add_option('-i', '--iterations', dest='iterations', default=1, help="The number of iterations to perform", type='int')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    for arg in args:
        cg = ftmc.CoarseGrainRNA(arg)
        cie = fbe.CylinderIntersectionEnergy()

        for i in range(options.iterations):
            sg = cg.random_subgraph()                                                                                        
            csg = ftmc.cg_from_sg(cg, sg)

            total_length = sum([len(list(cg.define_residue_num_iterator(d))) for d in sg])                                   
            cylinder_intersections = cie.calculate_intersection_coverages(csg)

            print(total_length, " ".join(map("{:.2f}".format, list(cylinder_intersections.values()))))

if __name__ == '__main__':
    main()

