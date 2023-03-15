#!/usr/bin/python

from __future__ import absolute_import
from __future__ import print_function
import sys, math
import warnings
from Bio.PDB import *

import tess.threedee.model.coarse_grain as ttmc
import borgy.graph.graph_pdb as cgg

def newest_output_graph(bg, chain):
    cgg.add_stem_information_from_pdb_chain(bg, chain)
    bg.add_bulge_coords_from_stems()
    cgg.add_loop_information_from_pdb_chain(bg, chain)


def main():
    if len(sys.argv) < 3:
        print("Usage: ./graph_to_simplified_pymol.py struct.graph pdb_file")
        sys.exit(1)

    pdb_name = sys.argv[2]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = PDBParser().get_structure('temp', pdb_name)
        chain = list(s.get_chains())[0]

    bg = ttmc.CoarseGrainRNA(sys.argv[1])

    cgg.add_stem_information_from_pdb_chain(bg, chain)
    bg.add_bulge_coords_from_stems()
    cgg.add_loop_information_from_pdb_chain(bg, chain)

    bg.dump()

if __name__=="__main__":
    main()
