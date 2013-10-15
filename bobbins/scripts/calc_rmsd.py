#!/usr/bin/python

import sys

from borgy.builder.rmsd import centered_rmsd
from borgy.graph.bulge_graph import BulgeGraph
import borgy.graph.graph_pdb as cgg

def main():
    if len(sys.argv) < 3:
        print >>sys.stderr, "Usage: ./calc_rmsd temp1.comp temp2.comp"
        print >>sys.stderr, "Calculate the rmsd between the two structures."
        sys.exit(1)

    bg1 = BulgeGraph(sys.argv[1])
    bg2 = BulgeGraph(sys.argv[2])

    #print "rmsd:", centered_rmsd(bg1.get_centers(), bg2.get_centers())
    print cgg.bg_rmsd(bg1, bg2)


if __name__ == '__main__':
    main()

