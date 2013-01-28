#!/usr/bin/python

import sys, os
import warnings
import numpy as np
import corgy.builder.stats as cbs
from optparse import OptionParser
import random
import corgy.builder.config as cbc
import corgy.builder.models as cbm
import corgy.builder.reconstructor as rtor
import corgy.utilities.debug as cud
import corgy.visual.pymol as cvp

import Bio.PDB as bpdb
import Bio.PDB.Chain as bpdbc

def main():
    usage = './align_stems.py [stem_length]'
    usage += 'Do diagnostics on the stem model'
    parser = OptionParser()

    parser.add_option('-i', '--iterations', dest='iterations', default=1, help="The number of times to repeat the alignment", type='int')
    parser.add_option('-l', '--length', dest='length', default=2, help="The length of the stem", type='int')
    parser.add_option('-o', '--output-pdb', dest='output_pdb', default=False, help="Output the structures to pdb files", action='store_true')

    (options, args) = parser.parse_args()

    if len(args) < 0:
        parser.print_help()
        sys.exit(1)

    stem_length = options.length
    if len(args) == 1:
        stem_length = int(args[0])

    sss = cbs.get_stem_stats()

    rmsds = []

    for i in range(options.iterations):
        stem_def = random.choice(sss[stem_length])
        filename = '%s_%s.pdb' % (stem_def.pdb_name, "_".join(map(str, stem_def.define)))
        pdb_file = os.path.join(cbc.Configuration.stem_fragment_dir, filename)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            chain = list(bpdb.PDBParser().get_structure('temp', pdb_file).get_chains())[0]

        m = cbm.define_to_stem_model(chain, stem_def.define)

        new_chain = bpdbc.Chain(' ')
        new_stem_def = random.choice(sss[stem_length])
        #cud.pv('filename')
        cbm.reconstruct_stem_core(new_stem_def, stem_def.define, new_chain, dict(), m)

        if options.output_pdb:
            rtor.output_chain(chain, 'out1.pdb')
            rtor.output_chain(new_chain, 'out3.pdb')

        unsuperimposed_rmsd = rtor.pdb_rmsd(chain, new_chain, backbone=False, superimpose=False)
        superimposed_rmsd = rtor.pdb_rmsd(chain, new_chain, backbone=False, superimpose=True)

        rmsds += [[superimposed_rmsd[1], unsuperimposed_rmsd[1]]]

        #cud.pv('(superimposed_rmsd, unsuperimposed_rmsd)')

        if options.output_pdb:
            rtor.output_chain(new_chain, 'out2.pdb')
            pp = cvp.PymolPrinter()
            (p,n) = m.mids
            pp.add_stem_like_core(m.mids, m.twists, stem_length+1, '')
            pp.stem_atoms(m.mids, m.twists, stem_length+1)
            pp.dump_pymol_file('ss')

    means = np.mean(np.array(rmsds), axis=0) 
    print stem_length, " ".join(map(str, means)), means[1] / means[0]

if __name__ == '__main__':
    main()

