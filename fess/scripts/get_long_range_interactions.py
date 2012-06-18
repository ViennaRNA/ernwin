#!/usr/bin/python

from Bio.PDB import *
import os, sys

from corgy.graph.bulge_graph import BulgeGraph
from corgy.utilities.mcannotate import iterate_over_interactions, get_interacting_base_pairs

def main():
    if len(sys.argv) < 3:
        print "Usage: ./get-long-range-interactions.py file.graph file.mcannotate"
        print 
        print "Print which graph elements interact with which others"
        exit

    bg = BulgeGraph(sys.argv[1])
    mcannotate_file = open(sys.argv[2], 'r')

    for line in iterate_over_interactions(mcannotate_file.readlines()):
        (from_chain, from_base, to_chain, to_base) =  get_interacting_base_pairs(line)

        node1 = bg.get_node_from_residue_num(from_base)
        node2 = bg.get_node_from_residue_num(to_base)

        if abs(from_base - to_base) > 1 and node1 != node2 and not bg.has_connection(node1, node2): # and not (node1[0] == 's' or node2[0] == 's'):
            bg.longrange[node1].add(node2)
            bg.longrange[node2].add(node1)

    bg.dump()

if __name__ == "__main__":
    main()
