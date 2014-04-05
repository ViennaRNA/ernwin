#!/usr/bin/python

import sys
import warnings

import Bio.PDB as bpdb
import borgy.builder.config as cbc
import borgy.builder.rmsd as cbr
import tess.threedee.model.coarse_grain as ttmc
import borgy.graph.graph_pdb as cgg
import borgy.utilities.debug as cud

from optparse import OptionParser

def main():
    usage = """
    ./virtual_atom_real_atom_distances.py temp.comp temp.pdb

    Calculate the distances between the virtual atoms and the real
    atoms.
    """
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        c = list(bpdb.PDBParser().get_structure('test', args[1]).get_chains())[0]

    if len(args) < 2:
        parser.print_help()
        sys.exit(1)

    bg = ttmc.CoarseGrainRNA(args[0])
    vposs = []
    rposs = []
    stem_length = 0
    for s in bg.stems():
        stem_length += bg.stem_length(s)
        cgg.add_virtual_residues(bg, s)
        for i in range(bg.stem_length(s)):
            for strand in range(2):
                va = cgg.virtual_residue_atoms(bg, s, i, strand)

                r = bg.stem_side_vres_to_resn(s, strand, i)
                for a, p in va.items():
                    try:
                        vec1 = p
                        vec2 = c[r][a].get_vector().get_array()

                        vposs += [vec1]
                        rposs += [vec2]
                    except KeyError as ke:
                        continue

    print cbc.Configuration.mids_method, stem_length, cbr.centered_rmsd(vposs, rposs)

if __name__ == '__main__':
    main()
