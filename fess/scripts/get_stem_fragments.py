#!/usr/bin/python

import corgy.graph.bulge_graph as cgb
import corgy.builder.config as cbc
import corgy.builder.stats as cbs

from Bio.PDB import PDBIO, PDBParser, Select

import os

import sys
from optparse import OptionParser

def output_stem_fragment(define, s, out_file):
    class StemSelect(Select):
        def accept_residue(self, residue):
            i = residue.id[1]

            if (define[0] <= i and i <= define[1]) or (define[2] <= i and i <= define[3]):
                return 1
            return 0

                    
    io=PDBIO()
    io.set_structure(s)
    io.save(out_file, StemSelect())

    print out_file

def main():
    usage = './get_stem_fragments.py [temp.comp]'
    parser = OptionParser()

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) == 1:
        bg = cgb.BulgeGraph(args[0])
        for st in bg.stems():
            filename = '%s_%s.pdb' % (bg.name, "_".join(map(str, bg.defines[st])))
            out_file = os.path.join(cbc.Configuration.stem_fragment_dir, filename)
            s = PDBParser().get_structure('t', os.path.join(cbc.Configuration.data_base_dir, "%s/prepare/temp.pdb" % (bg.name)))
            output_stem_fragment(bg.defines[st], s, out_file)
        sys.exit(0)

    stem_stats = cbs.get_stem_stats()

    structures = dict()
    prev_pdb_name = ''

    for sl in stem_stats.values():
        for ss in sl:
            filename = '%s_%s.pdb' % (ss.pdb_name, "_".join(map(str, ss.define)))
            out_file = os.path.join(cbc.Configuration.stem_fragment_dir, filename)

            if ss.pdb_name != prev_pdb_name:
                s = PDBParser().get_structure('t', os.path.join(cbc.Configuration.data_base_dir, "%s/prepare/temp.pdb" % (ss.pdb_name)))
                prev_pdb_name = ss.pdb_name

            output_stem_fragment(ss.define, s, out_file)


if __name__ == '__main__':
    main()

