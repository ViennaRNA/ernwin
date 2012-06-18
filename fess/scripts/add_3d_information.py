#!/usr/bin/python

import sys, math
from Bio.PDB import *

from corgy.graph.bulge_graph import BulgeGraph
from corgy.graph.graph_pdb import get_mids, get_bulge_centroid
from corgy.utilities.mcannotate import iterate_over_interactions, get_interacting_base_pairs
from corgy.utilities.vector import vec_distance

def output_stems(bg, chain):
    for d in bg.defines.keys():
        if d[0] == 's':
            mids = get_mids(chain, bg.defines[d])
            bg.coords[d] = (mids[0], mids[1])
            #print_segment(mids[0].get_array(), mids[1].get_array(), "green", 2.4, d)

def output_bulges(bg, chain):
    for d in bg.defines.keys():
        if d[0] != 's':
            if len(bg.edges[d]) == 2:
                edges = list(bg.edges[d])

                s1d = bg.defines[edges[0]]
                s2d = bg.defines[edges[1]]
                bd = bg.defines[d]
                
                (s1b, s1e) = bg.get_sides(edges[0], d)
                (s2b, s2e) = bg.get_sides(edges[1], d)

                mids1 = get_mids(chain, s1d)
                mids2 = get_mids(chain, s2d)

                #print >>sys.stderr, "mids1[s1b]:", mids1[s1b], mids2[s2b]
                if bg.weights[d] == 2:
                    bg.coords[d] = (mids1[s1b].get_array(), mids2[s2b].get_array())
                else:
                    bg.coords[d] = (mids1[s1b].get_array(), mids2[s2b]. get_array())


def output_loops(bg, chain):
    for d in bg.defines.keys():
        if d[0] != 's':
            if len(bg.edges[d]) == 1:
                s1 = list(bg.edges[d])[0]
                s1d = bg.defines[s1]
                bd = bg.defines[d]
                
                (s1b, s2b) = bg.get_sides(s1, d)

                mids = get_mids(chain, s1d)
                centroid = get_bulge_centroid(chain, bd)
                bg.coords[d] = (mids[s1b].get_array(), centroid.get_array())

def newest_output_graph(bg, chain):
    output_stems(bg, chain)
    output_bulges(bg, chain)
    output_loops(bg, chain)


def main():
    if len(sys.argv) < 3:
        print "Usage: ./graph_to_simplified_pymol.py struct.graph pdb_file"
        sys.exit(1)

    pdb_name = sys.argv[2]
    s = PDBParser().get_structure('temp', pdb_name)
    chain = list(s.get_chains())[0]

    bg = BulgeGraph(sys.argv[1])

    output_stems(bg, chain)
    output_bulges(bg, chain)
    output_loops(bg, chain)

    bg.dump()

if __name__=="__main__":
    main()
