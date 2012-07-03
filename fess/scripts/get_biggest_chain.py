#!/usr/bin/python

from Bio.PDB import *

import pdb
import sys

def main():
    if len(sys.argv) < 3:
        print "Usage: ./get_first_chain.py input_pdb output_pdb"
        sys.exit(1)

    s = PDBParser().get_structure('temp', sys.argv[1])

    chains = list(s.get_chains())
    biggest = 0
    biggest_len = 0

    for i in range(len(chains)):
        c = chains[i]
        num_residues = len(list(c.get_list()))

        if num_residues > biggest_len:
            biggest = i
            biggest_len = num_residues

    orig_chain = chains[biggest]

    class FirstChainSelect(Select):
        def accept_chain(self, chain):
            if chain.get_id() == orig_chain.get_id():
                return 1
            else:
                return 0

    io=PDBIO()
    io.set_structure(s)
    io.save(sys.argv[2], FirstChainSelect())

if __name__ == "__main__":
    main()
