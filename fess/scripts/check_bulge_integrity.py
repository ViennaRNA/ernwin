#!/usr/bin/python

from __future__ import absolute_import
from __future__ import print_function
from borgy.graph.bulge_graph import BulgeGraph
import math
import sys
from six.moves import range

def check_define_integrity(bg):
    for key in bg.defines.keys():
        d = bg.defines[key]

        for i in range(0, len(d), 2):
            #print >>sys.stderr, "key: %s i: %d d: %s" % (key, i, str(d))
            if d[i] > d[i+1]:
                print("PROBLEM: %s: d[%d] (%d) > d[%d] (%d)" % ( key, i, d[i], i+1, d[i+1] ), file=sys.stderr)
                sys.exit(1)
            #assert(d[i] <= d[i+1])

        if key[0] == 's':
            if d[1] - d[0] != d[3] - d[2]:
                print("PROBLEM: unequal stem lengths, s: %d d: %s" % (key, str(d)), file=sys.stderr)
            """
            if d[1] - d[0] == 0:
                print >> sys.stderr, "PROBLEM: stem length 1: %s" % (key)
            """

def check_connection_integrity(bg):
    nodes = set()

    for key in bg.defines.keys():
        nodes.add(key)

    for key in bg.edges.keys():
        for edge in bg.edges[key]:
            if key in nodes:
                nodes.remove(key)
            if edge in nodes:
                nodes.remove(edge)

    if len(nodes) > 0:
        print("PROBLEM: unconnected nodes:", file=sys.stderr)
        print("nodes", nodes, file=sys.stderr)
        sys.exit(1)


def main():
    if len(sys.argv) < 2:
        print("Usage: ./check_bulge_integrity.py struct.graph")
        print()
        print("Make sure the bulge graph satisfies the following constraints:")
        print("ds[0] < ds[1]")
        print("ds[2] < ds[3]")
        sys.exit(1)

    bg = BulgeGraph(sys.argv[1])
    check_define_integrity(bg)
    check_connection_integrity(bg)

if __name__=="__main__":
    main()
