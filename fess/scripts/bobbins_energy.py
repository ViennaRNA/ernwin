#!/usr/bin/python

import sys
from borgy.builder.energy import LongRangeDistanceEnergy
from borgy.graph.bulge_graph import BulgeGraph

def main():
    if len(sys.argv) < 2:
        print >>sys.stderr, "Usage: ./bobbins_energy.py temp.comp"
        sys.exit(1)

    bg = BulgeGraph(sys.argv[1])
    lrde = LongRangeDistanceEnergy()
    
    print lrde.naive_energy(bg), lrde.eval_energy(bg)

if __name__ == '__main__':
    main()

