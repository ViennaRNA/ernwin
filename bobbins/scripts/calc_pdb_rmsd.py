#!/usr/bin/python

import sys
import warnings

from optparse import OptionParser

import Bio.PDB as bpdb

import borgy.utilities.pdb as cup


def main():
    parser = OptionParser()

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < 2:
        print >>sys.stderr, "Usage: ./calc_pdb_rmsd.py pdb_file1 pdb_file2"
        sys.exit(1)

    fn1 = args[0]
    fn2 = args[1]

    rmsd = cup.pdb_file_rmsd(fn1, fn2)
    print rmsd[0], rmsd[1]

if __name__ == '__main__':
    main()

