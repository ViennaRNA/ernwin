#!/usr/bin/python

from corgy.builder.config import Configuration
from corgy.builder.stats import get_stem_stats

from Bio.PDB import PDBIO, PDBParser, Select

import os

stem_stats = get_stem_stats()

structures = dict()
prev_pdb_name = ''

for sl in stem_stats.values():
    for ss in sl:
        filename = '%s_%s.pdb' % (ss.pdb_name, "_".join(map(str, ss.define)))
        out_file = os.path.join(Configuration.stem_fragment_dir, filename)

        if ss.pdb_name != prev_pdb_name:
            s = PDBParser().get_structure('t', os.path.join(Configuration.data_base_dir, "%s/prepare/temp.pdb" % (ss.pdb_name)))
            prev_pdb_name = ss.pdb_name

        class StemSelect(Select):
            def accept_residue(self, residue):
                i = residue.id[1]

                if (ss.define[0] <= i and i <= ss.define[1]) or (ss.define[2] <= i and i <= ss.define[3]):
                    return 1
                return 0

                        
        io=PDBIO()
        io.set_structure(s)
        io.save(out_file, StemSelect())

        print out_file
