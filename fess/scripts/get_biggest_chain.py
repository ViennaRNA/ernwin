#!/usr/bin/python

from Bio.PDB import *

import pdb
import sys, warnings

def main():
    if len(sys.argv) < 3:
        print "Usage: ./get_first_chain.py input_pdb output_pdb"
        sys.exit(1)


    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if sys.argv[1] == '-':
            s = PDBParser().get_structure('temp', sys.stdin)
        else:
            s = PDBParser().get_structure('temp', sys.argv[1])

    chains = list(s.get_chains())
    biggest = 0
    biggest_len = 0

    for i in range(len(chains)):
        c = chains[i]
        res_list = list(c.get_list())
        
        #print >> sys.stderr, res_list[0].resname
        rna = False

        # Only count RNA residues
        num_residues = 0
        for res in res_list:
            if (res.resname.strip() == 'A' or
                res.resname.strip() == 'C' or
                res.resname.strip() == 'G' or
                res.resname.strip() == 'U'):
                num_residues += 1
                break

        if num_residues > biggest_len:
            biggest = i
            biggest_len = num_residues

    #print >>sys.stderr, "biggest", biggest
    orig_chain = chains[biggest]

    class FirstChainSelect(Select):
        def accept_chain(self, chain):
            if chain.get_id() == orig_chain.get_id():
                return 1
            else:
                return 0

    io=PDBIO()
    io.set_structure(s)

    if sys.argv[2] == '-':
        io.save(sys.stdout, FirstChainSelect())
    else:
        io.save(sys.argv[2], FirstChainSelect())

if __name__ == "__main__":
    main()
