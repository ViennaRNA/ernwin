import unittest
import corgy.aux.CPDB.src.examples.BarnacleCPDB as barn
import Bio.PDB as bpdb

import corgy.utilities.pdb as cup
import corgy.utilities.debug as cud

class TestPDB(unittest.TestCase):
    def test_num_noncovalent_clashes(self):
        loop_seq = 'AACCGGUUAAACCCGGGUUU'

        model = barn.BarnacleCPDB(loop_seq, 2.)
        model.sample()
        chain_loop = list(model.structure.get_chains())[0]

        cud.pv('cup.num_noncovalent_clashes(chain_loop)')
