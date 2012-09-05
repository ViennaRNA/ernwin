import unittest

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
from corgy.builder.ccd import get_closer_rotation_matrix

from corgy.builder.reconstructor import reconstruct_stems
from corgy.builder.reconstructor import output_chain
from corgy.builder.reconstructor import define_to_stem_model
from corgy.builder.reconstructor import splice_stem

from numpy import allclose, array, pi, dot, cross
from numpy.linalg import inv
from copy import deepcopy
from random import uniform, randint, random
from scipy.stats import norm, poisson

from Bio.PDB import PDBParser, PDBIO, Selection
from Bio.PDB.Chain import Chain
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure

import os

def get_alignment_vectors_rosetta(ress, r1, r2):
    #print "rosetta r1, r2:", r1, r2
    return( ress[r1]['O4*'].get_vector().get_array(),
            ress[r1]['C1*'].get_vector().get_array(),
            ress[r1]['C2*'].get_vector().get_array())

def get_alignment_vectors_barnacle(ress, r1, r2):
    #print "barnacle r1, r2:", r1, r2
    return( ress[r1]['O4\''].get_vector().get_array(),
            ress[r1]['C1\''].get_vector().get_array(),
            ress[r1]['C2\''].get_vector().get_array())

def get_measurement_vectors_rosetta(ress, r1, r2):
    return( ress[r2]['O4*'].get_vector().get_array(), 
            ress[r2]['C1*'].get_vector().get_array(),
            ress[r2]['C2*'].get_vector().get_array())

def get_measurement_vectors_barnacle(ress, r1, r2):
    return( ress[r2]['O4\''].get_vector().get_array(), 
            ress[r2]['C1\''].get_vector().get_array(),
            ress[r2]['C2\''].get_vector().get_array() )


def get_measurement_rmsd(chain1, chain2, handles):
    v1_m = get_measurement_vectors_rosetta(chain1, handles[0], handles[1])
    v2_m = get_measurement_vectors_barnacle(chain2, handles[2], handles[3])

    rmsd = vector_set_rmsd(v1_m, v2_m)

    return rmsd


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

            #print "stem.mids:", stem.mids
            #print "bg.coords:", sm.bg.coords[stem_name]
            #print

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
                stem = define_to_stem_model(chain, bg.defines[d])
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
                    #print start, end
                    #print stem.mids
                    #print
                    if ((allclose(stem.mids[0], start, atol=0.1) and allclose(stem.mids[1], end, atol=0.1)) or 
                        (allclose(stem.mids[1], start, atol=0.1) and allclose(stem.mids[0], end, atol=0.1))):
                            found += 1
                            break

        #print "FOUND:", found, "num_stems:", num_stems
        self.assertEquals(found, num_stems)

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

    def test_output_chain(self):
        bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))
        sm = SpatialModel(bg)

        sm.traverse_and_build()
        chain = reconstruct_stems(sm)

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

    def test_barnacle(self):
        sys.path.append('aux/Barnacle')
        from Barnacle import Barnacle

        model = Barnacle('ACGU')
        model.sample()

    def test_get_handles(self): 
        bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))
        sm = SpatialModel(bg)

        sm.traverse_and_build()
        chain = reconstruct_stems(sm)

        #fr = bg.get_flanking_region('b15', 0)
        (a,b,i1,i2) = bg.get_flanking_handles('b15')

        get_alignment_vectors_rosetta(chain, a, b)

    def test_barnacle_handles(self):
        sys.path.append('aux/Barnacle')
        from Barnacle import Barnacle

        bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))
        seq = bg.get_flanking_sequence('b15')
        (a,b,i1,i2) = bg.get_flanking_handles('b15')

        model = Barnacle(seq)
        model.sample()
        s = model.structure

        ress = list(s.get_residues())
        get_alignment_vectors_barnacle(ress, i1, i2)

    def test_handle_rmsd(self):
        sys.path.append('aux/Barnacle')
        from Barnacle import Barnacle

        bg = BulgeGraph(os.path.join(Configuration.test_input_dir, "1gid/graph", "temp.comp"))
        seq = bg.get_flanking_sequence('b15')
        (a,b,i1,i2) = bg.get_flanking_handles('b15')

        print "seq:", seq
        print "a:", a, "b:", "i1:", i1, "i2:", i2
        model = Barnacle(seq)
        model.sample()
        s = model.structure

        ress = list(s.get_residues())
        v1 = get_alignment_vectors_barnacle(ress, i1, i2)

        sm = SpatialModel(bg)
        sm.traverse_and_build()
        chain = reconstruct_stems(sm)

        print "a:", a, "b:", b
        v2 = get_alignment_vectors_rosetta(chain, a, b)
        prev_r = 1000.
        min_r = 1000.
        iterations = 10

        for i in range(iterations):
            sample_len = poisson.rvs(2)
            #print "sample_len:", sample_len

            if sample_len == 0:
                sample_len = 1

            j = randint(0, len(seq)-sample_len)
            m1 = j
            m2 = m1 + sample_len

            model.sample(start=m1, end=m2)

            ress = list(model.structure.get_residues())
            v1 = get_alignment_vectors_barnacle(ress, i1, i2)
            #r = centered_rmsd(v1, v2)
            r = self.quality_measurement(chain, list(model.structure.get_chains())[0], (a,b,i1,i2))
            print "r:", r

            #print "p:", p, "prev_p:", prev_p, "transition_p:", transition_p
            multiplier = .001 ** (1 / float(iterations))
            #print "multiplier:", multiplier
            temperature = 1. * (multiplier) ** i

            if r > prev_r:
                factor = -(r - prev_r) / ( temperature)
                
                p = exp (factor)
                #print "r-prev_r:", r-prev_r, "temperature:", temperature, "factor:", factor, "p:", p
                
                if random() > p:
                    model.undo()
                    continue

            print "r:", r, model.get_log_likelihood()
            prev_r = r

            if r < min_r:
                self.align_stems_and_barnacle(chain, list(model.structure.get_chains())[0], (a,b,i1,i2))
                min_r = r
            #sys.exit(1)

            self.assertGreater(centered_rmsd, 0)
            model.save_structure(os.path.join(Configuration.test_output_dir, 'best.pdb'))

    def perturb_c3_o3(self, chain, res_start, res_end, fixed):
        axis1 = chain[res_start]['C3\''].get_vector().get_array() - chain[res_start]['O3\''].get_vector().get_array()
        point1 = chain[res_start]['O3\''].get_vector().get_array()

        moving = get_measurement_vectors_barnacle(chain, res_start, res_end)

        rot_mat = get_closer_rotation_matrix(axis1, point1, moving, fixed)
        last_res = chain.get_list()[-1].id[1]

        for i in range(res_start+1, last_res+1):
            atoms = Selection.unfold_entities(chain[i], 'A')
            for atom in atoms:
                atom.transform(identity_matrix, -point1)
                atom.transform(rot_mat.transpose(), null_array)
                atom.transform(identity_matrix, point1)

    def perturb_p_o5(self, chain, res_start, res_end, fixed):
        axis1 = chain[res_start]['P'].get_vector().get_array() - chain[res_start]['O5\''].get_vector().get_array()
        point1 = chain[res_start]['O5\''].get_vector().get_array()

        moving = get_measurement_vectors_barnacle(chain, res_start, res_end)

        rot_mat = get_closer_rotation_matrix(axis1, point1, moving, fixed)
        last_res = chain.get_list()[-1].id[1]

        for i in range(res_start, last_res+1):
            atoms = Selection.unfold_entities(chain[i], 'A')
            for atom in atoms:
                atom.transform(identity_matrix, -point1)
                atom.transform(rot_mat.transpose(), null_array)
                atom.transform(identity_matrix, point1)

    def perturb_o5_c5(self, chain, res_start, res_end, fixed):
        axis1 = chain[res_start]['O5\''].get_vector().get_array() - chain[res_start]['C5\''].get_vector().get_array()
        point1 = chain[res_start]['C5\''].get_vector().get_array()

        moving = get_measurement_vectors_barnacle(chain, res_start, res_end)

        rot_mat = get_closer_rotation_matrix(axis1, point1, moving, fixed)
        last_res = chain.get_list()[-1].id[1]

        for i in range(res_start, last_res+1):
            atoms = Selection.unfold_entities(chain[i], 'A')
            for atom in atoms:
                if i == res_start and atom.name == 'P':
                    continue

                atom.transform(identity_matrix, -point1)
                atom.transform(rot_mat.transpose(), null_array)
                atom.transform(identity_matrix, point1)

    def perturb_c5_c4(self, chain, res_start, res_end, fixed):
        axis1 = chain[res_start]['C5\''].get_vector().get_array() - chain[res_start]['C4\''].get_vector().get_array()
        point1 = chain[res_start]['C4\''].get_vector().get_array()
        
        #print "axis1:", axis1
        #print "point1:", point1

        moving = get_measurement_vectors_barnacle(chain, res_start, res_end)

        #print "moving:", moving
        #print "fixed:", fixed

        rot_mat = get_closer_rotation_matrix(axis1, point1, moving, fixed)
        last_res = chain.get_list()[-1].id[1]

        for i in range(res_start, last_res+1):
            atoms = Selection.unfold_entities(chain[i], 'A')
            for atom in atoms:
                if i == res_start and (atom.name == 'P' or atom.name == 'O5\''):
                    continue

                atom.transform(identity_matrix, -point1)
                atom.transform(rot_mat.transpose(), null_array)
                atom.transform(identity_matrix, point1)

    def quality_measurement(self, chain_stems, chain_barnacle, handles):
        v1 = get_alignment_vectors_rosetta(chain_stems, handles[0], handles[1])
        v2 = get_alignment_vectors_barnacle(chain_barnacle, handles[2], handles[3])

        v1_m = get_measurement_vectors_rosetta(chain_stems, handles[0], handles[1])
        v2_m = get_measurement_vectors_barnacle(chain_barnacle, handles[2], handles[3])

        v1_centroid = get_vector_centroid(v1)
        v2_centroid = get_vector_centroid(v2)

        v1_t = v1 - v1_centroid
        v2_t = v2 - v2_centroid

        sup = optimal_superposition(v2_t, v1_t)

        v2_mt = v2_m - v2_centroid
        v2_rmt = dot(sup.transpose(), v2_mt.transpose()).transpose()

        v1_rmt = v1_m - v1_centroid

        rmsd = vector_set_rmsd(v1_rmt, v2_rmt)
        #print "rmsd1:", rmsd

        for atom in Selection.unfold_entities(chain_barnacle, 'A'):
            atom.transform(identity_matrix, -v2_centroid)
            atom.transform(sup, null_array)
            atom.transform(identity_matrix, v1_centroid)

        v1_m = get_measurement_vectors_rosetta(chain_stems, handles[0], handles[1])
        v2_m = get_measurement_vectors_barnacle(chain_barnacle, handles[2], handles[3])

        rmsd2 = get_measurement_rmsd(chain_stems, chain_barnacle, handles)
        print "rmsd2:", rmsd2

        for i in range(5):
            self.perturb_c3_o3(chain_barnacle, handles[2], handles[3], v1_m)

            for j in range(handles[2] + 1, handles[2] + 2):
                self.perturb_p_o5(chain_barnacle, j, handles[3], v1_m)
                self.perturb_o5_c5(chain_barnacle, j, handles[3], v1_m)
                self.perturb_c5_c4(chain_barnacle,  j, handles[3], v1_m)

            self.perturb_p_o5(chain_barnacle, handles[3], handles[3], v1_m)
            self.perturb_o5_c5(chain_barnacle, handles[3], handles[3], v1_m)
            self.perturb_c5_c4(chain_barnacle, handles[3], handles[3], v1_m)
            #print "i:", i, "rmsd9:", rmsd9

            '''
            self.assertGreater(rmsd2, rmsd3)
            self.assertGreater(rmsd3, rmsd4)
            self.assertGreater(rmsd4, rmsd5)
            self.assertGreater(rmsd5, rmsd6)
            '''
        rmsd9 = get_measurement_rmsd(chain_stems, chain_barnacle, handles)
        #print "axis1:", axis1, "point1:", point1
        return rmsd9

    def align_stems_and_barnacle(self, chain_stems, chain_barnacle, handles):
        v1 = get_alignment_vectors_rosetta(chain_stems, handles[0], handles[1])
        v2 = get_alignment_vectors_barnacle(chain_barnacle, handles[2], handles[3])

        v1_m = get_measurement_vectors_rosetta(chain_stems, handles[0], handles[1])
        v2_m = get_measurement_vectors_barnacle(chain_barnacle, handles[2], handles[3])

        v1_centroid = get_vector_centroid(v1)
        v2_centroid = get_vector_centroid(v2)

        v1_t = v1 - v1_centroid
        v2_t = v2 - v2_centroid
        sup = optimal_superposition(v1_t, v2_t)

        chain1 = Chain(' ')
        chain2 = Chain(' ')
        chain3 = Chain(' ')

        '''
        for res in chain_stems:
            chain1.add(res)
        '''
        counter = 0
        new_chain_barnacle = Chain(' ')

        for i in range(handles[2], handles[3]+1):
            id1 = chain_barnacle[i].id
            chain_barnacle[i].id = (id1[0], handles[0] + counter, id1[2])
            new_chain_barnacle.add(chain_barnacle[i])
            counter += 1

        chain_barnacle = new_chain_barnacle

        chain1.add(chain_stems[handles[0]])
        chain1.add(chain_stems[handles[1]])
        
        chain1.add(chain_stems[handles[0]-1])
        chain1.add(chain_stems[handles[1]+1])

        chain1.add(chain_stems[handles[0]-2])
        chain1.add(chain_stems[handles[1]+2])

        for i in range(handles[0], handles[1]+1):
            chain2.add(chain_barnacle[i])

        for atom in Selection.unfold_entities(chain1, 'A'):
            atom.transform(identity_matrix, -v1_centroid)
            atom.transform(sup, null_array)
            atom.transform(identity_matrix, v2_centroid)

        chain3.add(chain1[handles[0]-1])
        chain3.add(chain1[handles[1]+1])

        chain3.add(chain1[handles[0]-2])
        chain3.add(chain1[handles[1]+2])

        for i in range(handles[0], handles[1]+1):
            chain3.add(chain2[i])

        s1 = Structure(' ')
        m1 = Model(' ')
        m1.add(chain1)
        s1.add(m1)

        s2 = Structure(' ')
        m2 = Model(' ')
        m2.add(chain2)
        s2.add(m2)

        s3 = Structure(' ')
        m3 = Model(' ')
        m3.add(chain3)
        s3.add(m3)

        io = PDBIO()
        io.set_structure(s1)
        io.save(os.path.join(Configuration.test_output_dir, "s1.pdb"))

        io.set_structure(s2)
        io.save(os.path.join(Configuration.test_output_dir, 's2.pdb'))

        io.set_structure(s3)
        io.save(os.path.join(Configuration.test_output_dir, 's3.pdb'))


