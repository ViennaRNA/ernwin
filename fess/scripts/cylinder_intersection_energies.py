#!/usr/bin/python

import sys
from optparse import OptionParser

import forgi.threedee.model.coarse_grain as ftmc
import fess.builder.energy as cbe

def main():
    usage = './cylinder_intersection_energies.py bulge_graph'
    parser = OptionParser()

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    for arg in args:
        bg = ftmc.CoarseGrainRNA(arg)
        cie = cbe.CylinderIntersectionEnergy()
        
        in_cyl_fractions = cie.calculate_intersection_coverages(bg)
        
        for (key,val) in in_cyl_fractions.items():
            print key, val

if __name__ == '__main__':
    main()

