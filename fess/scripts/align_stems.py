#!/usr/bin/python

from __future__ import absolute_import
from __future__ import print_function
import sys, os
import warnings
import numpy as np
import borgy.builder.stats as cbs
from optparse import OptionParser
import random
import borgy.builder.config as cbc
import borgy.builder.models as cbm
import borgy.builder.stats as cbs
import borgy.builder.reconstructor as rtor
import borgy.graph.graph_pdb as cgg
import borgy.utilities.debug as cud
import borgy.utilities.pdb as cup
import borgy.visual.pymol as cvp

import Bio.PDB as bpdb
import Bio.PDB.Chain as bpdbc
from six.moves import map
from six.moves import range

def stem_def_from_filename(filename):
    parts = filename.split('.')[0].split('_')

    stem_def = cbs.StemStat()
    stem_def.pdb_name = parts[0]
    stem_def.define = [int(parts[1]), int(parts[2]), int(parts[3]), int(parts[4])]
    stem_def.bp_length = stem_def.define[1] - stem_def.define[0]

    return stem_def

def main():
    usage = './align_stems.py [stem_length]'
    usage += 'Do diagnostics on the stem model'
    parser = OptionParser()

    parser.add_option('-i', '--iterations', dest='iterations', default=1, help="The number of times to repeat the alignment", type='int')
    parser.add_option('-l', '--length', dest='length', default=2, help="The length of the stem", type='int')
    parser.add_option('-o', '--output-pdb', dest='output_pdb', default=False, help="Output the structures to pdb files", action='store_true')
    parser.add_option('-f', '--from', dest='from_file', default=None, help='Specify a file to align from. Invalidates the -l option.', type='str')
    parser.add_option('-t', '--to', dest='to_file', default=None, help='Specify a file to align to. Invalidates the -l option.', type='str')
    parser.add_option('-m', '--method', dest='method', default='e', help='Specify which method to use for the helix fitting. e = estimate (original, least accurate method), a = align (better, more accurate method), t = template (best, most accurate method)')
    parser.add_option('-a', '--average-twist', dest='use_average_method', default=False, action='store_true', help='Use the average of the two twists to align the stems.')

    (options, args) = parser.parse_args()

    if len(args) < 0:
        parser.print_help()
        sys.exit(1)

    stem_length = options.length
    if len(args) == 1:
        stem_length = int(args[0])

    if options.from_file == None or options.to_file == None:
        sss = cbs.get_stem_stats(os.path.join(cbc.Configuration.base_dir, 'fess/stats/temp.1jj2.stats'))

    rmsds = []

    for i in range(options.iterations):
        if options.from_file != None:
            filename = options.from_file
            stem_def = stem_def_from_filename(filename)
        else:
            stem_def = random.choice(sss[stem_length])
            filename = '%s_%s.pdb' % (stem_def.pdb_name, "_".join(map(str, stem_def.define)))

        pdb_file = os.path.join(cbc.Configuration.stem_fragment_dir, filename)

        # Extract the PDB coordinates of the original chain
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            try:
                chain = list(bpdb.PDBParser().get_structure('temp', pdb_file).get_chains())[0]
                chain = cbm.extract_stem_from_chain(chain, stem_def)
            except IOError as ie:
                cud.pv('ie')

        # Convert the chain into a stem model
        # This is where the method for fitting a helix is applied
        #m = cbm.define_to_stem_model(chain, stem_def.define)
        stem = cbm.StemModel(name=stem_def.define)
        define = stem_def.define
        mids = cgg.get_mids(chain, define, options.method)

        stem.mids = tuple([m.get_array() for m in mids])
        stem.twists = cgg.get_twists(chain, define)
        m = stem

        # Create a new chain by aligning the stem from the sampled define
        # to the model created from the original stem
        new_chain = bpdbc.Chain(' ')
        try:
            if options.to_file != None:
                new_stem_def = stem_def_from_filename(options.to_file)
            else:
                new_stem_def = random.choice(sss[stem_def.bp_length])

            cbm.reconstruct_stem_core(new_stem_def, stem_def.define, new_chain, dict(), m, options.use_average_method)
        except IOError as ie:
            cud.pv('ie')

        if options.output_pdb:
            rtor.output_chain(chain, 'out1.pdb')
            rtor.output_chain(new_chain, 'out3.pdb')

        unsuperimposed_rmsd = cup.pdb_rmsd(chain, new_chain, sidechains=False, superimpose=False)
        superimposed_rmsd = cup.pdb_rmsd(chain, new_chain, sidechains=False, superimpose=True, apply_sup=True)
        rmsds += [[superimposed_rmsd[1], unsuperimposed_rmsd[1]]]

        #cud.pv('(superimposed_rmsd, unsuperimposed_rmsd)')

        if options.output_pdb:
            rtor.output_chain(new_chain, 'out2.pdb')
            pp = cvp.PymolPrinter()
            (p,n) = m.mids
            pp.add_stem_like_core(m.mids, m.twists, stem_def.bp_length+1, '')
            pp.stem_atoms(m.mids, m.twists, stem_def.bp_length+1)
            pp.dump_pymol_file('ss')

        print(stem_length, superimposed_rmsd[1], unsuperimposed_rmsd[1], unsuperimposed_rmsd[1] / superimposed_rmsd[1])

    #means = np.mean(np.array(rmsds), axis=0) 
    #print stem_length, " ".join(map(str, means)), means[1] / means[0]

if __name__ == '__main__':
    main()

