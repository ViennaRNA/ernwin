#!/usr/bin/python

from __future__ import absolute_import
from __future__ import print_function
import sys
from optparse import OptionParser

import borgy.utilities.pdb as cup
import Bio.PDB as bp
from six.moves import map

def main():
    usage = './noncovalent_distances.py pdb_file'
    usage += 'Print the distances between all atoms which are'
    usage += 'closer than 5 angstroms and not bonded covalently'

    parser = OptionParser()

    parser.add_option('-c', '--cutoff', dest='cutoff', default=3., help="The maximum distance for considering two atoms", type='float')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    s = bp.PDBParser().get_structure('temp', args[0])
    distances = cup.noncovalent_distances(list(s.get_chains())[0], options.cutoff)
    print("\n".join(map(str,distances)))

if __name__ == '__main__':
    main()

