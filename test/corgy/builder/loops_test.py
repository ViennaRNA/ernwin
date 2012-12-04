import unittest
import os

import corgy.builder.config as cbc
import corgy.builder.loops as cbl
import corgy.builder.models as cbm
import corgy.builder.reconstructor as rtor
import corgy.graph.bulge_graph as cgb

class TestLoops(unittest.TestCase):
    def test_build_loop(self):
        bg = cgb.BulgeGraph(os.path.join(cbc.Configuration.test_input_dir, "1y26/graph", "temp.comp"))
        sm = cbm.SpatialModel(bg)
        sm.sample_native_stems()
        sm.create_native_stem_models()

        #sm.traverse_and_build()
        chain = rtor.reconstruct_stems(sm)
        ld = 'b3'
        bg = sm.bg
        side = 0
        seq = bg.get_flanking_sequence(ld, side)
        (a,b,i1,i2) = bg.get_flanking_handles(ld, side)

        best_loop_chain = cbl.build_loop(chain, seq, (a,b,i1,i2), bg.length, 15, consider_contacts=True)

        cbl.add_loop_chain(chain, best_loop_chain, (a,b,i1,i2), bg.length)
        rtor.output_chain(chain, os.path.join(cbc.Configuration.test_output_dir, 'r1.pdb'))
