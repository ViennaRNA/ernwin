#!/usr/bin/python

import sys

import pickle 
import traceback

from sys import stderr

from optparse import OptionParser
from corgy.builder.energy import LongRangeInteractionCount, CombinedEnergy
from corgy.builder.energy import JunctionClosureEnergy, SkewNormalInteractionEnergy
from corgy.builder.models import SpatialModel

import corgy.builder.energy as cbe
import corgy.builder.config as conf

import os, pdb

from corgy.graph.bulge_graph import BulgeGraph
import numpy as np

def main():
    parser = OptionParser()
    np.seterr(all='ignore')

    energies = []
    energies += [SkewNormalInteractionEnergy()]
    energies += [LongRangeInteractionCount()]
    energies += [JunctionClosureEnergy()]

    #uncal_energies = [cbe.StemClashEnergy()]
    uncal_energies = [cbe.StemVirtualResClashEnergy()]

    #energies += [LongRangeDistanceEnergy()]

    parser.add_option('-i', '--iterations', dest='iterations', default=1000, help="The number of calibration steps.", type='int')
    parser.add_option('-e', '--energy-dir', dest='energy_dir', default=os.path.join(conf.Configuration.base_dir, 'bobbins/energy'), help="The root energy saving directory", type='string')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        print >>sys.stderr, "Usage: ./calibrate_energy temp.comp"
        print >>sys.stderr, "Create an energy function, calibrate it to the structure in temp.comp"
        print >>sys.stderr, "Then dump it to a file."
        sys.exit(1)


    ce = CombinedEnergy(energies, uncal_energies)

    bg = BulgeGraph(args[0])

    ce.calibrate(SpatialModel(bg), iterations=options.iterations, bg_energy=None, output_dir=options.energy_dir)

if __name__ == '__main__':
    try:
        main()
    except:
        tb = traceback.format_exc()
        print >>stderr, tb
        #pdb.post_mortem()
        raise

