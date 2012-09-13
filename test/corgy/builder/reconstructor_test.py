import unittest

import aux.Barnacle as barn

import corgy.builder.reconstructor as rc
import pdb, sys

from operator import attrgetter
from corgy.visual.pymol import PymolPrinter
from math import sqrt, exp

from pprint import pprint
from corgy.builder.config import Configuration

from corgy.graph.bulge_graph import BulgeGraph
from corgy.builder.models import StemModel

from corgy.graph.graph_pdb import get_mids, get_twists
from corgy.graph.graph_pdb import get_stem_orientation_parameters
from corgy.utilities.vector import rotation_matrix, magnitude
from corgy.utilities.vector import x_array, null_array, identity_matrix
from corgy.utilities.vector import get_vector_centroid

from corgy.builder.models import SpatialModel
from corgy.builder.rmsd import centered_rmsd, optimal_superposition

import corgy.builder.reconstructor as rtor

from numpy import allclose, array, pi, dot, cross
from numpy.linalg import inv
from copy import deepcopy
from random import uniform, random

from Bio.PDB import PDBParser, PDBIO, Selection
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure

import Bio.PDB as bpdb
import corgy.utilities.vector as cuv

import os

def get_random_orientation():
    '''
    Return a random tuple (u, v, t) such that
    0 <= u <= pi
    -pi <= v <= pi
    -pi <= t <= pi
    '''

    return (uniform(0, pi), uniform(-pi, pi), uniform(-pi, pi))

def get_random_translation():
    '''
    Return a random translation.
    '''

    return array([uniform(-10, 10), uniform(-10, 10), uniform(-10, 10)])

class TestReconstructor(unittest.TestCase):
    def setUp(self):
        self.bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "2b3j/graph", "temp.comp"))
        s = PDBParser().get_structure('temp', os.path.join(Configuration.test_input_dir, "2b3j/prepare", "temp.pdb"))

        self.chain = list(s.get_chains())[0]
        self.stem = rtor.define_to_stem_model(self.chain, self.bg.defines['s0'])
    
    def test_rotate_stem(self):
        stem1 = StemModel()

        stem1.mids = self.bg.coords['s0']
        stem1.twists = self.bg.twists['s0']

        stem2 = rtor.rotate_stem(stem1, get_random_orientation())

        self.assertFalse(allclose(stem1.twists[0], stem2.twists[0]))
        self.assertFalse(allclose(stem1.twists[1], stem2.twists[1]))

        self.assertTrue(allclose(stem1.mids[0], stem2.mids[0]))
        self.assertTrue(allclose(magnitude(stem1.mids[1] - stem1.mids[0]), magnitude(stem2.mids[1] - stem2.mids[0])))

    def test_define_to_stem_model(self):
        stem1 = StemModel()
        stem1.mids = self.bg.coords['s0']
        stem1.twists = self.bg.twists['s0']

        stem2 = rtor.define_to_stem_model(self.chain, self.bg.defines['s0'])

        self.assertTrue(stem1 == stem2)

    def test_rerotate_stem(self):
        stem1 = deepcopy(self.stem)

        orientation = get_random_orientation()
        stem2 = rtor.rotate_stem(stem1, orientation)

        # vector should not be rotated away from itself... duh!
        (r, u, v, t) = get_stem_orientation_parameters(stem1.vec(), stem1.twists[0], stem1.vec(), stem1.twists[0])
        self.assertTrue(allclose((u,v,t), (pi/2,0.,0.)))

        (r, u, v, t) = get_stem_orientation_parameters(stem1.vec(), stem1.twists[0], stem2.vec(), stem2.twists[0])

        # Figure out why exactly this works!!!
        orientation1 = (pi-u, -v, -t)
        rot_mat = rtor.get_stem_rotation_matrix(stem1, orientation1)
        
        stem3 = deepcopy(stem2)
        stem3.rotate(inv(rot_mat), offset=stem3.mids[0])

        (r, u, v, t) = get_stem_orientation_parameters(stem1.vec(), stem1.twists[0], stem3.vec(), stem3.twists[0])
        self.assertTrue(allclose((u,v,t),(pi/2, 0., 0.)))

    def test_splice_stem(self):
        define = self.bg.defines['s0']

        start1 = define[0]
        end1 = define[1]

        start2 = define[2]
        end2 = define[3]

        new_chain = rtor.splice_stem(self.chain, define)
        residues = Selection.unfold_entities(new_chain, 'R')

        # Make sure only the residues specified are in the chain
        for res in residues:
            self.assertTrue((res.id[1] >= start1 and res.id[1] <= end1) or (res.id[1] >= start2 and res.id[1] <= end2))

        # try to access each residue to see if they are really there
        # if not, the testing framework will catch the error and report it
        for i in xrange(start1, end1+1):
            res = new_chain[i]
        
        for i in xrange(start2, end2+1):
            res = new_chain[i]

    def test_rotate_atom_stem(self):
        chain = rtor.splice_stem(self.chain, self.bg.defines['s0'])

        stem1 = rtor.define_to_stem_model(chain, self.bg.defines['s0'])
        orientation = get_random_orientation()
        stem2 = rtor.rotate_stem(stem1, orientation)

        self.assertFalse(stem1 == stem2)

        (r, u, v, t) = get_stem_orientation_parameters(stem1.vec(), stem1.twists[0], stem2.vec(), stem2.twists[0])
        rot_mat = rtor.get_stem_rotation_matrix(stem1, (pi-u, -v, -t))
        rtor.rotate_chain(chain, inv(rot_mat), stem1.mids[0])

        stem3 = rtor.define_to_stem_model(chain, self.bg.defines['s0'])

        self.assertTrue(stem2 == stem3)

    def test_align_chain_to_stem(self):
        chain = rtor.splice_stem(self.chain, self.bg.defines['s0'])

        stem1 = rtor.define_to_stem_model(chain, self.bg.defines['s0'])
        orientation = get_random_orientation()
        translation = get_random_translation()

        stem2 = rtor.rotate_stem(stem1, orientation)
        #stem2.translate(translation)

        self.assertFalse(stem1 == stem2)
        rtor.align_chain_to_stem(chain, self.bg.defines['s0'], stem2)
        stem3 = rtor.define_to_stem_model(chain, self.bg.defines['s0'])

        self.assertTrue(stem2 == stem3)

    def test_align_chain_to_stem1(self):
        chain = rtor.splice_stem(self.chain, self.bg.defines['s0'])

        stem1 = rtor.define_to_stem_model(chain, self.bg.defines['s0'])
        orientation = get_random_orientation()
        translation = get_random_translation()

        stem2 = rtor.rotate_stem(stem1, orientation)
        stem2.translate(translation)

        self.assertFalse(stem1 == stem2)
        rtor.align_chain_to_stem(chain, self.bg.defines['s0'], stem2)
        stem3 = rtor.define_to_stem_model(chain, self.bg.defines['s0'])

        self.assertTrue(stem2 == stem3)

    def check_reconstructed_stems(self, sm, chain, stem_names):
        for stem_name in stem_names:
            stem_def = sm.stem_defs[stem_name]
            bg_stem_def = sm.bg.defines[stem_name]
        
            stem = rtor.define_to_stem_model(chain, bg_stem_def)

            self.assertEqual(stem, sm.stems[stem_name])
    
    def check_pymol_stems(self, bg, coarse_filename, pdb_filename):
        '''
        Check whether the file output by pymol_printer is consistent
        with the output pdb file.
        '''
        found = 0
        num_stems = 0

        chain = list(PDBParser().get_structure('t', pdb_filename).get_chains())[0]
        stems = []
        for d in bg.defines.keys():
            if d[0] == 's':
                stem = rtor.define_to_stem_model(chain, bg.defines[d])
                stems += [stem]

                num_stems += 1

        f = open(coarse_filename, 'r')
        cylinders = []

        for line in f:
            if line.find('CYLINDER') >= 0:
                parts = line.strip().split(', ')
                start = array(map(float, parts[1:4]))
                end = array(map(float, parts[4:7]))

                for stem in stems:
                    if ((allclose(stem.mids[0], start, atol=0.1) and allclose(stem.mids[1], end, atol=0.1)) or 
                        (allclose(stem.mids[1], start, atol=0.1) and allclose(stem.mids[0], end, atol=0.1))):
                            found += 1
                            break

        self.assertEquals(found, num_stems)

    def test_reconstruct_sampled_whole_model(self):
        '''
        Test the reconstruction of the stems of the SpatialModel.
        '''
        bgs = []
        #bgs += [self.bg]
        bgs += [BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))]

        for bg in bgs:
            sm = SpatialModel(bg)
            sm.traverse_and_build()
            chain = rtor.reconstruct_stems(sm)

            self.check_reconstructed_stems(sm, chain, sm.stem_defs.keys())

    def test_reconstruct_native_whole_model(self):
        bgs = []
        #bgs += [self.bg]
        bgs += [BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))]

        for bg in bgs:
            sm = SpatialModel(bg)
            sm.sample_native_stems()
            sm.create_native_stem_models()
            #sm.traverse_and_build()
            chain = rtor.reconstruct_stems(sm)

            self.check_reconstructed_stems(sm, chain, sm.stem_defs.keys())
    
    def test_twice_defined_stem(self):
        bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))
        sm = SpatialModel(bg)

        sm.traverse_and_build()
        sm.stem_defs['s2'].pdb_name = '1s72'
        sm.stem_defs['s2'].bp_length = 1
        sm.stem_defs['s2'].phys_length = 3.78678002926
        sm.stem_defs['s2'].twist_angle = 0.457289932542
        sm.stem_defs['s2'].define = [678, 679, 684, 685]

        sm.stem_defs['s1'].pdb_name = '1yjn'
        sm.stem_defs['s1'].bp_length = 1
        sm.stem_defs['s1'].phys_length = 3.77964149852
        sm.stem_defs['s1'].twist_angle = 0.455189553351
        sm.stem_defs['s1'].define = [678, 679, 684, 685]
        sm.traverse_and_build()

        chain = rtor.reconstruct_stems(sm)
        self.check_reconstructed_stems(sm, chain, sm.stem_defs.keys())

    def test_output_chain(self):
        bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))
        sm = SpatialModel(bg)

        sm.traverse_and_build()
        chain = rtor.reconstruct_stems(sm)

        self.check_reconstructed_stems(sm, chain, sm.stem_defs.keys())

        output_file = os.path.join(Configuration.test_output_dir, "output_chain")

        pymol_printer = PymolPrinter()
        pymol_printer.print_text = False
        pymol_printer.add_twists = True
        pymol_printer.add_longrange = False

        pymol_printer.coordinates_to_pymol(sm.bg)
        pymol_printer.chain_to_pymol(chain)
        pymol_printer.dump_pymol_file(output_file)

        chain = list(PDBParser().get_structure('t', output_file + ".pdb").get_chains())[0]

        self.check_reconstructed_stems(sm, chain, sm.stem_defs.keys())
        self.check_pymol_stems(bg, output_file + ".pym", output_file + ".pdb")

    def test_reconstruct_loops(self):
        bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))
        sm = SpatialModel(bg)
        sm.sample_native_stems()
        sm.create_native_stem_models()

        #sm.traverse_and_build()
        chain = rtor.reconstruct_stems(sm)
        rtor.reconstruct_loops(chain, sm)
        '''
        rtor.reconstruct_loop(chain, sm, 'b15')
        #rtor.reconstruct_loop(chain, sm, 'b1')
        rtor.reconstruct_loop(chain, sm, 'b11')
        rtor.reconstruct_loop(chain, sm, 'b18')
        #rtor.reconstruct_loop(chain, sm, 'b16')
        #rtor.reconstruct_loop(chain, sm, 'x2', 0)
        #rtor.reconstruct_loop(chain, sm, 'x2', 1)
        '''

        #self.check_reconstructed_stems(sm, chain, sm.stem_defs.keys())
        rtor.output_chain(chain, os.path.join(Configuration.test_output_dir, 'r1.pdb'))

    def test_reconstruct_loop(self):
        bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))
        sm = SpatialModel(bg)
        sm.sample_native_stems()
        sm.create_native_stem_models()

        #sm.traverse_and_build()
        chain = rtor.reconstruct_stems(sm)
        rtor.reconstruct_loop(chain, sm, 'b18')
        #rtor.reconstruct_loop(chain, sm, 'x4', side=1)
        '''
        #rtor.reconstruct_loop(chain, sm, 'b1')
        rtor.reconstruct_loop(chain, sm, 'b11')
        rtor.reconstruct_loop(chain, sm, 'b18')
        #rtor.reconstruct_loop(chain, sm, 'b16')
        #rtor.reconstruct_loop(chain, sm, 'x2', 0)
        #rtor.reconstruct_loop(chain, sm, 'x2', 1)
        '''

        #self.check_reconstructed_stems(sm, chain, sm.stem_defs.keys())
        rtor.output_chain(chain, os.path.join(Configuration.test_output_dir, 'r1.pdb'))

    def test_get_stem_coord_array(self):
        bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))
        ld = 'b15'
        seq = bg.get_flanking_sequence(ld)
        (a,b,i1,i2) = bg.get_flanking_handles(ld)
        model = barn.Barnacle(seq)
        model.sample()
        chain = list(model.structure.get_chains())[0]

        (coords, indeces) = rtor.get_atom_coord_array(chain, i1, i2)
        
        for i in range(i1, i2+1):
            res = chain[i]
            self.assertTrue(allclose(coords[indeces[res.id[1]]], res['P'].get_vector().get_array()))

    def test_close_fragment_loop(self):
        import aux.Barnacle as barn
        bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))
        sm = SpatialModel(bg)
        sm.sample_native_stems()
        sm.create_native_stem_models()

        ld = 'b15'
        seq = bg.get_flanking_sequence(ld)
        (a,b,i1,i2) = bg.get_flanking_handles(ld)

        model = barn.Barnacle(seq)
        model.sample()
        s = model.structure

        chain_stems = rtor.reconstruct_stems(sm) 
        chain_barnacle = list(model.structure.get_chains())[0]

        rtor.align_starts(chain_stems, chain_barnacle, (a,b,i1,i2))
        distances = rtor.get_adjacent_interatom_distances(chain_barnacle, i1, i2)
        (moving, indeces) = rtor.get_atom_coord_array(chain_barnacle, i1, i2)

        r, chain_loop = rtor.close_fragment_loop(chain_stems, chain_barnacle, (a,b,i1,i2), iterations=10)

        (moving1, indeces1) = rtor.get_atom_coord_array(chain_loop, i1, i2)

        # make sure that the loop closure actually did something
        self.assertFalse(allclose(moving, moving1))
        distances1 = rtor.get_adjacent_interatom_distances(chain_loop, i1, i2)
        self.assertTrue(allclose(distances, distances1))

    def test_barnacle(self):
        from aux.Barnacle import Barnacle

        model = Barnacle('ACGU')
        model.sample()

    def test_get_handles(self): 
        bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))
        sm = SpatialModel(bg)

        sm.traverse_and_build()
        chain = rtor.reconstruct_stems(sm)

        #fr = bg.get_flanking_region('b15', 0)
        (a,b,i1,i2) = bg.get_flanking_handles('b15')

        rtor.get_alignment_vectors(chain, a, b)

    def test_pdb_rmsd(self):
        s1 = bpdb.PDBParser().get_structure('t', os.path.join(Configuration.test_input_dir, '1gid/prepare/temp.pdb'))
        s2 = bpdb.PDBParser().get_structure('t', os.path.join(Configuration.test_input_dir, '1gid/prepare/temp_sampled.pdb'))

        sup = bpdb.Superimposer()

        rs1 = list(s1.get_residues())
        rs2 = list(s2.get_residues())

        c1 = list(s1.get_chains())[0]
        c2 = list(s2.get_chains())[0]

        a_5_names = ['P', 'O5*', 'C5*', 'C4*', 'O4*', 'O2*']
        a_3_names = ['C1*', 'C2*', 'C3*', 'O3*']

        a_names = dict()
        a_names['U'] = a_5_names + ['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6'] + a_3_names
        a_names['C'] = a_5_names + ['N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6'] + a_3_names

        a_names['A'] = a_5_names + ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N6', 'N7', 'C8', 'N9'] + a_3_names
        a_names['G'] = a_5_names + ['N1', 'C2', 'N2', 'N3', 'C4', 'C5', 'C6', 'O6', 'N7', 'C8', 'N9'] + a_3_names

        all_atoms1 = []
        all_atoms2 = []

        for i in range(1, len(list(s1.get_residues()))+1):
            #anames = a_5_names + a_names[c1[i].resname.strip()] + a_3_names
            anames = a_5_names + a_3_names

            atoms1 = [c1[i][a] for a in anames]
            atoms2 = [c2[i][a] for a in anames]

            if len(atoms1) != len(atoms2):
                print "different lengths"
                sys.exit(1)

            all_atoms1 += atoms1
            all_atoms2 += atoms2

        sup.set_atoms(all_atoms1, all_atoms2)

        print "sup.rms:", sup.rms


                

