#!/usr/bin/python

from __future__ import absolute_import
from __future__ import print_function
import sys
from borgy.graph.bulge_graph import BulgeGraph

def main():
    if len(sys.argv) < 2:
        print("Usage: ./bulge_to_fork_graph.py bulge.graph")

    bg = BulgeGraph(sys.argv[1])
    bg.collapse()
    bg.dump()

if __name__ == "__main__":
    main()
