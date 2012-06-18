#!/usr/bin/python

from Bio.PDB import *
import os, sys
from bulge_to_fork_graph import parse_graph, BulgeGraph
from graph_to_simplified_pymol import get_bulge_centroid, print_connection

def parse_chain_base(chain_base):
    #print >>sys.stderr,"chain_base:", chain_base
    if ord(chain_base[0]) >= ord('A') and ord(chain_base[0]) <= ord('Z'):
        chain = chain_base[0]
        base = int(chain_base[1:])
    else:
        if chain_base[0] == '\'':
            end_quote_idx = chain_base.find('\'',1)
            chain=chain_base[1:end_quote_idx]
            base=int(chain_base[end_quote_idx+1:])
        else:
            chain=' '
            base=int(chain_base)

    return (chain, base)

def parse_base_pair_id(base_pair_id):
    parts = base_pair_id.split('-')
    (from_chain, from_base) = parse_chain_base(parts[0].strip())
    (to_chain, to_base) = parse_chain_base(parts[1].strip())

    return (from_chain, from_base, to_chain, to_base)

def get_interacting_base_pairs(line):
    line_parts = line.split(' ')
    return parse_base_pair_id(line_parts[0])

def iterate_over_interactions(mcannotate_lines):
    base_pair_line = False
    for line in mcannotate_lines:
        if line.find("Base-pairs ---") == 0:
            base_pair_line = True
            continue 
        if line.find("Residue conformations") == 0:
            base_pair_line = False
            continue
        if base_pair_line:
            yield line.strip()

def main():
    if len(sys.argv) < 4:
        print "Usage: ./long-range-to-pymol.py file.graph file.mcannotate pdb_file"
        print 
        print "Print which graph elements interact with which others"
        exit

    bg = parse_graph(open(sys.argv[1], 'r'))
    mcannotate_file = open(sys.argv[2], 'r')

    pdb_name = sys.argv[3]
    s = PDBParser().get_structure('temp', pdb_name)
    c = list(s.get_chains())[0]

    print "from pymol.cgo import *"
    print "from pymol import cmd"
    print "from pymol.vfont import plain"
    print "obj = ["
    for line in iterate_over_interactions(mcannotate_file.readlines()):
        (from_chain, from_base, to_chain, to_base) =  get_interacting_base_pairs(line)

        node1 = bg.get_node_from_residue_num(from_base)
        node2 = bg.get_node_from_residue_num(to_base)

        if abs(from_base - to_base) > 1 and node1 != node2 and not bg.has_connection(node1, node2) and not (node1[0] == 's' or node2[0] == 's'):
            c1 = get_bulge_centroid(c, bg.defines[node1])
            c2 = get_bulge_centroid(c, bg.defines[node2])
            print_connection(c1, c2)

    print " ]"
    print "cmd.load_cgo(obj, 'ss')"
if __name__ == "__main__":
    main()
