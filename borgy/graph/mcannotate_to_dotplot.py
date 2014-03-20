#!/usr/bin/python

import sys
import copy

from get_long_range_interactions import iterate_over_interactions, iterate_over_residue_list, parse_chain_base

def get_dotplot(lines):
    """docstring for get_dotplot"""
    residues = []
    bps = DefaultDict('.')

    for line in iterate_over_residue_list(lines):
        parts = line.split(' ')
        residues += [parts[0]]

    for line in iterate_over_interactions(lines):
        parts = line.split(' ')
        bond_type = parts[3]
        if bond_type.find('Ww/Ww') >= 0 or bond_type.find('Ww/Ws') >= 0 or bond_type.find('Ws/Ww') >= 0:
            parts1 = parts[0].split('-')
            bps[parts1[0]] = '('
            bps[parts1[1]] = ')'

    out_str = ''
    for res in residues:
        out_str += bps[res]

    print out_str

    pass

def main():
    """docstring for main"""
    if len(sys.argv) < 2:
        print "Usage: ./mcannotate_to_dotplot.py"
        exit(1)

    lines = open(sys.argv[1], 'r').readlines()
    
    get_dotplot(lines)

    pass

if __name__ == '__main__':
    main()
