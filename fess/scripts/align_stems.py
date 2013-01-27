#!/usr/bin/python

import sys, os
import corgy.builder.stats as cbs
from optparse import OptionParser
import random
import corgy.builder.config as cbc
import corgy.builder.models as cbm
import corgy.builder.reconstructor as rtor
import corgy.utilities.debug as cud
import Bio.PDB as bpdb
import Bio.PDB.Chain as bpdbc

def main():
    usage = './align_stems.py [stem_length]'
    usage += 'Do diagnostics on the stem model'
    parser = OptionParser()

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < 0:
        parser.print_help()
        sys.exit(1)

    stem_length = 3
    if len(args) == 1:
        stem_length = int(args[0])

    sss = cbs.get_stem_stats()

    stem_def = random.choice(sss[stem_length])
    filename = '%s_%s.pdb' % (stem_def.pdb_name, "_".join(map(str, stem_def.define)))
    pdb_file = os.path.join(cbc.Configuration.stem_fragment_dir, filename)
    chain = list(bpdb.PDBParser().get_structure('temp', pdb_file).get_chains())[0]
    m = cbm.define_to_stem_model(chain, stem_def.define)

    new_chain = bpdbc.Chain(' ')
    new_stem_def = random.choice(sss[stem_length])
    cud.pv('filename')
    cbm.reconstruct_stem_core(new_stem_def, stem_def.define, new_chain, dict(), m)

    rtor.output_chain(chain, 'out1.pdb')
    rtor.output_chain(new_chain, 'out3.pdb')
    cud.pv('rtor.pdb_rmsd(chain, new_chain, backbone=False, superimpose=False)')
    cud.pv('rtor.pdb_rmsd(chain, new_chain, backbone=False, superimpose=True)')
    rtor.output_chain(new_chain, 'out2.pdb')


if __name__ == '__main__':
    main()

