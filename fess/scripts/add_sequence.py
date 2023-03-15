#!/usr/bin/python

from __future__ import absolute_import
from __future__ import print_function
import sys
from optparse import OptionParser

from borgy.graph.bulge_graph import BulgeGraph

def main():
    parser = OptionParser()

    (options, args) = parser.parse_args()

    if len(args) < 2:
        print("Usage: ./add_sequence.py struct.graph seq_file", file=sys.stderr)
        sys.exit(1)

    bg = BulgeGraph(args[0])

    f = open(args[1])
    lines = f.readlines()
    bg.seq = lines[0].strip()
    f.close()

    bg.dump()

if __name__ == '__main__':
    main()

