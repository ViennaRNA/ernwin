#!/usr/bin/python

import sys
import copy

from corgy.utilities.mcannotate import iterate_over_interactions, iterate_over_residue_list, parse_chain_base

# From:
# http://code.activestate.com/recipes/389639/
class DefaultDict(dict):
    """Dictionary with a default value for unknown keys."""
    def __init__(self, default):
        self.default = default

    def __getitem__(self, key):
        if key in self: 
            return self.get(key)
        else:
            ## Need copy in case self.default is something like []
            return self.setdefault(key, copy.deepcopy(self.default))

    def __copy__(self):
        copy = DefaultDict(self.default)
        copy.update(self)
        return copy

def get_dotplot(lines):
    """docstring for get_dotplot"""
    residues = []
    residue_types = []
    bps = DefaultDict(-1)

    for line in iterate_over_residue_list(lines):
        parts = line.split(' ')
        residues += [parts[0]]
        residue_types += [parts[2]]

    paired = set()
    for line in iterate_over_interactions(lines):
        parts = line.split(' ')
        bond_type = parts[3]
        if bond_type.find('Ww/Ww') >= 0 or bond_type.find('Ww/Ws') >= 0 or bond_type.find('Ws/Ww') >= 0:
        #if bond_type.find('Ww/Ww') >= 0:
            parts1 = parts[0].split('-')
            #print line

            if parts1[0] in paired or parts1[1] in paired:
                print >>sys.stderr, "paired:", parts1[0], parts1[1]
                continue

            paired.add(parts1[0])
            paired.add(parts1[1])

            bps[parts1[0]] = residues.index(parts1[1])
            bps[parts1[1]] = residues.index(parts1[0])

    for i in range(len(residue_types)):
        print i+1, residue_types[i], bps[residues[i]]+1

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
