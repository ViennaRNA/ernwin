#!/usr/bin/python

from __future__ import absolute_import
from __future__ import print_function
import collections as c
import itertools as it
import pdb

import forgi.threedee.model.coarse_grain as ttmc
import forgi.threedee.model.stats as cbs
import fess.builder.config as fbc
import forgi.utilities.debug as cud

from Bio.PDB import PDBIO, PDBParser, Select

import os

import sys
from optparse import OptionParser
from six.moves import map
from six.moves import range

def output_stem_fragment(define, s, out_file):
    class StemSelect(Select):
        def accept_residue(self, residue):
            i = residue.id[1]

            for j in range(0, len(define), 2):
                if define[j]-3 <= i and i <= define[j+1]+3:
                    return 1
            return 0

            '''
            if (define[0] <= i and i <= define[1]) or (define[2] <= i and i <= define[3]):
                return 1
            return 0
            '''

                    
    io=PDBIO()
    io.set_structure(s)
    io.save(out_file, StemSelect())

def main():
    usage = './get_stem_fragments.py [temp.comp]'
    parser = OptionParser()

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')
    parser.add_option('-o', '--output-dir', dest='output_dir', default=fbc.Configuration.stem_fragment_dir, 
                      help='The directory in which to output all of the fragments', type='str')

    (options, args) = parser.parse_args()

    if not os.path.exists(options.output_dir):
        os.makedirs(options.output_dir)

    if len(args) == 1:
        bg = ttmc.CoarseGrainRNA(args[0])

        for st in bg.elements():
            filename = '%s_%s.pdb' % (bg.name, "_".join(map(str, bg.defines[st])))
            out_file = os.path.join(options.output_dir, filename)
            s = PDBParser().get_structure('t', os.path.join(fbc.Configuration.data_base_dir, "%s/temp.pdb" % (bg.name)))
            output_stem_fragment(bg.defines[st], s, out_file)
        sys.exit(0)

    #stats = [cbs.get_angle_stats(), cbs.get_loop_stats()]
    #stem_stats = cbs.get_stem_stats()

    structures = dict()
    prev_pdb_name = ''

    for l in it.chain(cbs.get_angle_stats().values(), cbs.get_stem_stats().values(),
                       cbs.get_loop_stats().values()):
        for ss in l:
            filename = '%s_%s.pdb' % (ss.pdb_name, "_".join(map(str, ss.define)))
            out_file = os.path.join(options.output_dir, filename)

            if ss.pdb_name != prev_pdb_name:
                cud.pv('ss.define, fbc.Configuration.data_base_dir, ss.pdb_name')
                s = PDBParser().get_structure('t', os.path.join(fbc.Configuration.data_base_dir, "%s/temp.pdb" % (ss.pdb_name)))
                prev_pdb_name = ss.pdb_name

            print(out_file, ss.define)
            output_stem_fragment(ss.define, s, out_file)


if __name__ == '__main__':
    main()

