#!/usr/bin/python

import sys
import warnings

from optparse import OptionParser

import Bio.PDB as bpdb

import corgy.utilities.pdb as cup

def main():
    parser = OptionParser()

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < 2:
        print >>sys.stderr, "Usage: ./calc_pdb_rmsd.py pdb_file1 pdb_file2"
        sys.exit(1)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        s1= bpdb.PDBParser().get_structure('t', args[0])
        s2= bpdb.PDBParser().get_structure('t', args[1])

    c1 = list(s1.get_chains())[0]
    c2 = list(s2.get_chains())[0]

    rmsd = cup.pdb_rmsd(c1, c2)
    print rmsd[0], rmsd[1]

if __name__ == '__main__':
    main()

