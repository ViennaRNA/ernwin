#!/usr/bin/python

from Bio.PDB import *

import sys

def main():
    if len(sys.argv) < 2:
        print "Usage: ./separate_all_chains.py pdb_file"
        sys.exit(1)

    s = PDBParser().get_structure('temp', sys.argv[1])

    chains = list(s.get_chains())
    orig_chain = list(s.get_chains())[0]
    counter = 0

    for c in chains:
        class FirstChainSelect(Select):
            def accept_chain(self, chain):
                if chain.get_id() == c.get_id():
                    return 1
                else:
                    return 0

        io = PDBIO()
        io.set_structure(s)
        io.save(sys.argv[1]+"."+str(counter), FirstChainSelect())

        counter += 1

if __name__ == "__main__":
    main()
