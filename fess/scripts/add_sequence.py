#!/usr/bin/python

import sys
from optparse import OptionParser

from corgy.graph.bulge_graph import BulgeGraph

def main():
    parser = OptionParser()

    (options, args) = parser.parse_args()

    if len(args) < 2:
        print >> sys.stderr, "Usage: ./add_sequence.py struct.graph seq_file"
        sys.exit(1)

    bg = BulgeGraph(args[0])

    f = open(args[1])
    lines = f.readlines()
    bg.seq = lines[0].strip()
    f.close()

    bg.dump()

if __name__ == '__main__':
    main()

