import unittest

import corgy.builder.reconstructor as rc
from operator import attrgetter

from pprint import pprint
from corgy.builder.config import Configuration

from corgy.graph.bulge_graph import BulgeGraph
from corgy.builder.models import StemModel

from corgy.graph.graph_pdb import get_mids, get_twists
from corgy.graph.graph_pdb import get_stem_orientation_parameters
from corgy.utilities.vector import rotation_matrix, magnitude
from corgy.utilities.vector import x_array, null_array, identity_matrix

from corgy.builder.models import SpatialModel

from numpy import allclose, array, pi, dot, cross
from numpy.linalg import inv
from copy import deepcopy
from random import uniform


from Bio.PDB import PDBParser, PDBIO, Selection
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure

import os

def reconstruct_stems(sm):
    '''
    Reconstruct the stems around a Spatial Model.

    @param sm: Spatial Model
    '''
    #sm.traverse_and_build()
    new_chain = Chain(' ')

    for stem_name in sm.stem_defs.keys():
        stem_def = sm.stem_defs[stem_name]
        stem = sm.stems[stem_name]

        filename = '%s_%s.pdb' % (stem_def.pdb_name, "_".join(map(str, stem_def.define)))
        #print "stem_name:", stem_name, "stem_def:", stem_def
        pdb_file = os.path.join(Configuration.stem_fragment_dir, filename)

        chain = list(PDBParser().get_structure('temp', pdb_file).get_chains())[0]
        align_chain_to_stem(chain, stem_def.define, stem)

        #print "stem_def.define:", stem_def.define
        #print "sm.bg.defines[stem_name]:", sm.bg.defines[stem_name]

        #for e in chain.get_list():
        for i in range(stem_def.bp_length+1):
            #print "i:", i
            e = chain[stem_def.define[0] + i]
            e.id = (e.id[0], sm.bg.defines[stem_name][0] + i, e.id[2])
            #print "adding:", e.id
            new_chain.add(e)

            e = chain[stem_def.define[2] + i]
            e.id = (e.id[0], sm.bg.defines[stem_name][2] + i, e.id[2])
            #print "adding:", e.id
            new_chain.add(e)

    return new_chain

def define_to_stem_model(chain, define):
    '''
    Extract a StemModel from a Bio.PDB.Chain structure.

    The define is 4-tuple containing the start and end coordinates
    of the stem on each strand. 

    s1s s1e s2s s2e

    @param chain: The Bio.PDB.Chain representation of the chain
    @param define: The BulgeGraph define
    @return: A StemModel with the coordinates and orientation of the stem.
    '''
    stem = StemModel()

    mids = get_mids(chain, define)

    stem.mids = tuple([m.get_array() for m in mids])
    stem.twists = get_twists(chain, define)

    return stem

def splice_stem(chain, define):
    '''
    Extract just the defined stem from the chain and return it as
    a new chain.
    
    @param chain: A Bio.PDB.Chain containing the stem in define
    @param define: The BulgeGraph stem define
    '''
    start1 = define[0]
    end1 = define[1]

    start2 = define[2]
    end2 = define[3]

    new_chain = Chain(' ')

    for i in xrange(start1, end1+1):
        #new_chain.insert(i, chain[i])
        new_chain.add(chain[i])

    for i in xrange(start2, end2+1):
        new_chain.add(chain[i])

    '''
    m = Model(' ')
    s = Structure(' ')
    m.add(new_chain)
    s.add(m)

    io=PDBIO()
    io.set_structure(s)
    io.save('temp.pdb')
    '''

    return new_chain

def get_stem_rotation_matrix(stem, (u, v, t)):
    twist1 = stem.twists[0]

    # rotate around the stem axis to adjust the twist


    # rotate down from the twist axis
    comp1 = cross(stem.vec(), twist1)

    rot_mat1 = rotation_matrix(stem.vec(), t)
    rot_mat2 = rotation_matrix(twist1, u - pi/2)
    rot_mat3 = rotation_matrix(comp1, v)

    rot_mat4 = dot(rot_mat3, dot(rot_mat2, rot_mat1))

    return rot_mat4

def rotate_stem(stem, (u, v, t)):
    '''
    Rotate a particular stem.
    '''
    stem2 = deepcopy(stem)
    rot_mat4 = get_stem_rotation_matrix(stem, (u,v,t))
    stem2.rotate(rot_mat4, offset=stem.mids[0])

    return stem2

def rotate_chain(chain, rot_mat, offset):
    atoms = Selection.unfold_entities(chain, 'A')

    for atom in atoms:
        atom.transform(identity_matrix, -offset)
        atom.transform(rot_mat, null_array)
        atom.transform(identity_matrix, offset)

def translate_chain(chain, translation):
    atoms = Selection.unfold_entities(chain, 'A')

    for atom in atoms:
        atom.transform(identity_matrix, translation)


def align_chain_to_stem(chain, define, stem2):
    stem1 = define_to_stem_model(chain, define)

    (r, u, v, t) = get_stem_orientation_parameters(stem1.vec(), stem1.twists[0], stem2.vec(), stem2.twists[0])
    rot_mat = get_stem_rotation_matrix(stem1, (pi-u, -v, -t))
    rotate_chain(chain, inv(rot_mat), stem1.mids[0])
    translate_chain(chain, stem2.mids[0] - stem1.mids[0])


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
        self.stem = define_to_stem_model(self.chain, self.bg.defines['s0'])
    
    def test_rotate_stem(self):
        stem1 = StemModel()

        stem1.mids = self.bg.coords['s0']
        stem1.twists = self.bg.twists['s0']

        stem2 = rotate_stem(stem1, get_random_orientation())

        self.assertFalse(allclose(stem1.twists[0], stem2.twists[0]))
        self.assertFalse(allclose(stem1.twists[1], stem2.twists[1]))

        self.assertTrue(allclose(stem1.mids[0], stem2.mids[0]))
        self.assertTrue(allclose(magnitude(stem1.mids[1] - stem1.mids[0]), magnitude(stem2.mids[1] - stem2.mids[0])))

    def test_define_to_stem_model(self):
        stem1 = StemModel()
        stem1.mids = self.bg.coords['s0']
        stem1.twists = self.bg.twists['s0']

        stem2 = define_to_stem_model(self.chain, self.bg.defines['s0'])

        self.assertTrue(stem1 == stem2)

    def test_rerotate_stem(self):
        stem1 = deepcopy(self.stem)

        orientation = get_random_orientation()
        stem2 = rotate_stem(stem1, orientation)

        # vector should not be rotated away from itself... duh!
        (r, u, v, t) = get_stem_orientation_parameters(stem1.vec(), stem1.twists[0], stem1.vec(), stem1.twists[0])
        self.assertTrue(allclose((u,v,t), (pi/2,0.,0.)))

        (r, u, v, t) = get_stem_orientation_parameters(stem1.vec(), stem1.twists[0], stem2.vec(), stem2.twists[0])

        # Figure out why exactly this works!!!
        orientation1 = (pi-u, -v, -t)
        rot_mat = get_stem_rotation_matrix(stem1, orientation1)
        
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

        new_chain = splice_stem(self.chain, define)
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
        chain = splice_stem(self.chain, self.bg.defines['s0'])

        stem1 = define_to_stem_model(chain, self.bg.defines['s0'])
        orientation = get_random_orientation()
        stem2 = rotate_stem(stem1, orientation)

        self.assertFalse(stem1 == stem2)

        (r, u, v, t) = get_stem_orientation_parameters(stem1.vec(), stem1.twists[0], stem2.vec(), stem2.twists[0])
        rot_mat = get_stem_rotation_matrix(stem1, (pi-u, -v, -t))
        rotate_chain(chain, inv(rot_mat), stem1.mids[0])

        stem3 = define_to_stem_model(chain, self.bg.defines['s0'])

        self.assertTrue(stem2 == stem3)

    def test_align_chain_to_stem(self):
        chain = splice_stem(self.chain, self.bg.defines['s0'])

        stem1 = define_to_stem_model(chain, self.bg.defines['s0'])
        orientation = get_random_orientation()
        translation = get_random_translation()

        stem2 = rotate_stem(stem1, orientation)
        #stem2.translate(translation)

        self.assertFalse(stem1 == stem2)
        align_chain_to_stem(chain, self.bg.defines['s0'], stem2)
        stem3 = define_to_stem_model(chain, self.bg.defines['s0'])

        self.assertTrue(stem2 == stem3)

    def test_align_chain_to_stem1(self):
        chain = splice_stem(self.chain, self.bg.defines['s0'])

        stem1 = define_to_stem_model(chain, self.bg.defines['s0'])
        orientation = get_random_orientation()
        translation = get_random_translation()

        stem2 = rotate_stem(stem1, orientation)
        stem2.translate(translation)

        self.assertFalse(stem1 == stem2)
        align_chain_to_stem(chain, self.bg.defines['s0'], stem2)
        stem3 = define_to_stem_model(chain, self.bg.defines['s0'])

        self.assertTrue(stem2 == stem3)

    def check_reconstructed_stems(self, sm, chain, stem_names):
        for stem_name in stem_names:
            stem_def = sm.stem_defs[stem_name]
            bg_stem_def = sm.bg.defines[stem_name]
        
            stem = define_to_stem_model(chain, bg_stem_def)

            print "stem_def:", stem_def
            print stem
            print sm.stems[stem_name]
            print

            self.assertEqual(stem, sm.stems[stem_name])

    '''
    def test_reconstruct_one_stem(self):
        sm = SpatialModel(self.bg)
        sm.traverse_and_build()
        stem_def = sm.stem_defs['s0']

        filename = '%s_%s.pdb' % (stem_def.pdb_name, "_".join(map(str, stem_def.define)))
        pdb_file = os.path.join(Configuration.stem_fragment_dir, filename)

        chain = list(PDBParser().get_structure('temp', pdb_file).get_chains())[0]
        align_chain_to_stem(chain, stem_def.define, sm.stems['s0'])
        stem = define_to_stem_model(chain, sm.stem_defs['s0'].define)

        #self.assertTrue(stem == sm.stems['s0'])

        self.check_reconstructed_stems(sm, chain, ['s0'])
    '''

    def test_reconstruct_whole_model(self):
        '''
        Test the reconstruction of the stems of the SpatialModel.
        '''
        bgs = []
        #bgs += [self.bg]
        bgs += [BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))]

        for bg in bgs:
            sm = SpatialModel(bg)
            sm.traverse_and_build()
            chain = reconstruct_stems(sm)

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

        chain = reconstruct_stems(sm)
        self.check_reconstructed_stems(sm, chain, sm.stem_defs.keys())
