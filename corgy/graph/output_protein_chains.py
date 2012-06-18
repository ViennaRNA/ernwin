#!/usr/bin/python

import sys
from Bio.PDB import *

def main():
    if len(sys.argv) < 3:
        print "Usage: ./output-protein-chain.py pdb_file output_file"
        print
        print "Output a PDB file consisting of only protein chains. In other words, remove RNA chains."
        sys.exit(1)

    struct = PDBParser().get_structure('temp', sys.argv[1])

    class ProteinSelect(Select):
        def accept_chain(self, chain):
            for res in chain:
                if res.resname.strip() in ['A', 'U', 'C', 'G', 'rA', 'rU', 'rC', 'rG']:
                    return 0
            return 1

    io = PDBIO()
    io.set_structure(struct)
    io.save(sys.argv[2], ProteinSelect())

    
if __name__ == "__main__":
    main()
