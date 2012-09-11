import corgy.builder.models as models
import corgy.builder.rmsd as brmsd

import corgy.graph.graph_pdb as gpdb
import corgy.utilities.vector as cuv

import corgy.builder.ccd as cbc

import aux.Barnacle as barn

import Bio.PDB as bpdb
import Bio.PDB.Chain as bpdbc
import Bio.PDB.Model as bpdbm
import Bio.PDB.Structure as bpdbs

import scipy.stats as ss

from scipy.stats import norm, poisson

import os, math, sys
import corgy.builder.config as conf
import copy
import random as rand

import numpy as np

def get_alignment_vectors_rosetta(ress, r1, r2):
    return( ress[r1]['O4*'].get_vector().get_array(),
            ress[r1]['C1*'].get_vector().get_array(),
            ress[r1]['C2*'].get_vector().get_array())

def get_alignment_vectors_barnacle(ress, r1, r2):
    return( ress[r1]['O4*'].get_vector().get_array(),
            ress[r1]['C1*'].get_vector().get_array(),
            ress[r1]['C2*'].get_vector().get_array())

def get_measurement_vectors_rosetta(ress, r1, r2):
    return( ress[r2]['O4*'].get_vector().get_array(), 
            ress[r2]['C1*'].get_vector().get_array(),
            ress[r2]['C2*'].get_vector().get_array())

def get_measurement_vectors_barnacle(ress, r1, r2):
    return( ress[r2]['O4*'].get_vector().get_array(), 
            ress[r2]['C1*'].get_vector().get_array(),
            ress[r2]['C2*'].get_vector().get_array() )


def get_measurement_rmsd(chain1, chain2, handles):
    v1_m = get_measurement_vectors_rosetta(chain1, handles[0], handles[1])
    v2_m = get_measurement_vectors_barnacle(chain2, handles[2], handles[3])

    rmsd = cuv.vector_set_rmsd(v1_m, v2_m)

    return rmsd
def get_stem_rotation_matrix(stem, (u, v, t)):
    twist1 = stem.twists[0]

    # rotate around the stem axis to adjust the twist

    # rotate down from the twist axis
    comp1 = np.cross(stem.vec(), twist1)

    rot_mat1 = cuv.rotation_matrix(stem.vec(), t)
    rot_mat2 = cuv.rotation_matrix(twist1, u - math.pi/2)
    rot_mat3 = cuv.rotation_matrix(comp1, v)

    rot_mat4 = np.dot(rot_mat3, np.dot(rot_mat2, rot_mat1))

    return rot_mat4

def rotate_stem(stem, (u, v, t)):
    '''
    Rotate a particular stem.
    '''
    stem2 = copy.deepcopy(stem)
    rot_mat4 = get_stem_rotation_matrix(stem, (u,v,t))
    stem2.rotate(rot_mat4, offset=stem.mids[0])

    return stem2

def rotate_chain(chain, rot_mat, offset):
    '''
    Move according to rot_mat for the position of offset.

    @param chain: A Bio.PDB.Chain instance.
    @param rot_mat: A left_multiplying rotation_matrix.
    @param offset: The position from which to do the rotation.
    '''

    atoms = bpdb.Selection.unfold_entities(chain, 'A')

    for atom in atoms:
        #atom.transform(np.eye(3,3), -offset)
        atom.coord -= offset
        atom.transform(rot_mat, offset)

def translate_chain(chain, translation):
    '''
    Translate all of the atoms in a chain by a certain amount.

    @param chain: A Bio.PDB.Chain instance to be translated.
    @translation: A vector indicating the direction of the translation.
    '''
    atoms = bpdb.Selection.unfold_entities(chain, 'A')

    for atom in atoms:
        atom.transform(np.eye(3,3), translation)

def align_chain_to_stem(chain, define, stem2):
    stem1 = define_to_stem_model(chain, define)

    (r, u, v, t) = gpdb.get_stem_orientation_parameters(stem1.vec(), stem1.twists[0], stem2.vec(), stem2.twists[0])
    rot_mat = get_stem_rotation_matrix(stem1, (math.pi-u, -v, -t))
    rotate_chain(chain, np.linalg.inv(rot_mat), stem1.mids[0])
    translate_chain(chain, stem2.mids[0] - stem1.mids[0])

def reconstruct_stems(sm):
    '''
    Reconstruct the stems around a Spatial Model.

    @param sm: Spatial Model
    '''
    #sm.traverse_and_build()
    new_chain = bpdbc.Chain(' ')

    for stem_name in sm.stem_defs.keys():
        stem_def = sm.stem_defs[stem_name]
        stem = sm.stems[stem_name]

        filename = '%s_%s.pdb' % (stem_def.pdb_name, "_".join(map(str, stem_def.define)))
        #print "stem_name:", stem_name, "stem_def:", stem_def
        pdb_file = os.path.join(conf.Configuration.stem_fragment_dir, filename)

        chain = list(bpdb.PDBParser().get_structure('temp', pdb_file).get_chains())[0]
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


def output_chain(chain, filename):
    '''
    Dump a chain to an output file.

    @param chain: The Bio.PDB.Chain to dump.
    @param filename: The place to dump it.
    '''
    m = bpdbm.Model(' ')
    s = bpdbs.Structure(' ')

    m.add(chain)
    s.add(m)

    io = bpdb.PDBIO()
    io.set_structure(s)
    io.save(filename)

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
    stem = models.StemModel()

    mids = gpdb.get_mids(chain, define)

    stem.mids = tuple([m.get_array() for m in mids])
    stem.twists = gpdb.get_twists(chain, define)

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

    new_chain = bpdbc.Chain(' ')

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

def add_residue_to_rosetta_chain(chain, residue):
    '''
    Add a residue and rename all of it's atoms to the Rosetta convention.

    C1' -> C1*

    @param chain: The chain to add to
    @param residue: The residue to be added
    '''
    removed_atoms = []

    for atom in residue.get_list():
        removed_atoms += [atom]
        residue.detach_child(atom.id)

        atom.name = atom.name.replace('\'', '*')
        atom.id = atom.name

    for atom in removed_atoms:
        residue.add(atom)

    chain.add(residue)

def add_loop_chain(chain, loop_chain, handles, length):
    '''
    Add all of the residues in loop_chain to chain.

    @param chain: The target chain to which the residues will be added.
    @param loop_chain: The source of the loop residues.
    @param handles: The indeces of the adjacent stem regions as well as the indeces into the loop
        chain which define which section is actually the loop and which is the additional linker
        region.
    '''
    # detach the residues of the helix which are adjacent to the loop
    #r1_id = chain[handles[0]].id
    #chain.detach_child(r1_id)
    #replace them with the residues of the loop
    #loop_chain[handles[2]].id = r1_id
    #add_residue_to_rosetta_chain(chain, loop_chain[handles[2]])

    if handles[1] != length:
        r2_id = chain[handles[1]].id
        chain.detach_child(r2_id)
        loop_chain[handles[3]].id = r2_id
        add_residue_to_rosetta_chain(chain, loop_chain[handles[3]])

    counter = 1
    for i in range(handles[2]+1, handles[3]):
        id1 = loop_chain[i].id
        loop_chain[i].id = (id1[0], handles[0] + counter, id1[2])
        add_residue_to_rosetta_chain(chain, loop_chain[i])
        counter += 1


def perturb_c3_o3(chain, res_start, res_end, fixed):
    axis1 = chain[res_start]['C3\''].get_vector().get_array() - chain[res_start]['O3\''].get_vector().get_array()
    point1 = chain[res_start]['O3\''].get_vector().get_array()

    moving = get_measurement_vectors_barnacle(chain, res_start, res_end)

    rot_mat = cbc.get_closer_rotation_matrix(axis1, point1, moving, fixed)
    last_res = chain.get_list()[-1].id[1]

    for i in range(res_start+1, last_res+1):
        atoms = bpdb.Selection.unfold_entities(chain[i], 'A')
        for atom in atoms:
            atom.transform(np.eye(3,3), -point1)
            atom.transform(rot_mat.transpose(), point1)

def perturb_p_o5(chain, res_start, res_end, fixed):
    axis1 = chain[res_start]['P'].get_vector().get_array() - chain[res_start]['O5\''].get_vector().get_array()
    point1 = chain[res_start]['O5\''].get_vector().get_array()

    moving = get_measurement_vectors_barnacle(chain, res_start, res_end)

    rot_mat = cbc.get_closer_rotation_matrix(axis1, point1, moving, fixed)
    last_res = chain.get_list()[-1].id[1]

    for i in range(res_start, last_res+1):
        atoms = bpdb.Selection.unfold_entities(chain[i], 'A')
        for atom in atoms:
            atom.transform(np.eye(3,3), -point1)
            atom.transform(rot_mat.transpose(), point1)

def perturb_o5_c5(chain, res_start, res_end, fixed):
    axis1 = chain[res_start]['O5\''].get_vector().get_array() - chain[res_start]['C5\''].get_vector().get_array()
    point1 = chain[res_start]['C5\''].get_vector().get_array()

    moving = get_measurement_vectors_barnacle(chain, res_start, res_end)

    rot_mat = cbc.get_closer_rotation_matrix(axis1, point1, moving, fixed)
    last_res = chain.get_list()[-1].id[1]

    for i in range(res_start, last_res+1):
        atoms = bpdb.Selection.unfold_entities(chain[i], 'A')
        for atom in atoms:
            if i == res_start and atom.name == 'P':
                continue

            atom.transform(np.eye(3,3), -point1)
            atom.transform(rot_mat.transpose(), point1)

def perturb_c5_c4(chain, res_start, res_end, fixed):
    axis1 = chain[res_start]['C5\''].get_vector().get_array() - chain[res_start]['C4\''].get_vector().get_array()
    point1 = chain[res_start]['C4\''].get_vector().get_array()
    
    moving = get_measurement_vectors_barnacle(chain, res_start, res_end)

    rot_mat = cbc.get_closer_rotation_matrix(axis1, point1, moving, fixed)
    last_res = chain.get_list()[-1].id[1]

    for i in range(res_start, last_res+1):
        atoms = bpdb.Selection.unfold_entities(chain[i], 'A')
        for atom in atoms:
            if i == res_start and (atom.name == 'P' or atom.name == 'O5\''):
                continue

            atom.transform(np.eye(3,3), -point1)
            atom.transform(rot_mat.transpose(), point1)

a_5_names = ['P', 'O5*', 'C5*', 'C4*', 'O4*', 'C1*', 'C2*'] 
a_3_names = ['C3*', 'O3*']

a_names = dict()
a_names['U'] = a_5_names + ['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6'] + a_3_names
a_names['C'] = a_5_names + ['N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6'] + a_3_names

a_names['A'] = a_5_names + ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N6', 'N7', 'C8', 'N9'] + a_3_names
a_names['G'] = a_5_names + ['N1', 'C2', 'N2', 'N3', 'C4', 'C5', 'C6', 'O6', 'N7', 'C8', 'N9'] + a_3_names

def get_atom_coord_array(chain, start_res, end_res):
    '''
    Return an array of the coordinates of all the atoms in the chain, 
    arranged in the following order:

    P, O5', C5', C4', O4', C3', C1', C2', base_atoms, O3'

    @param chain: The chain from which to get the coordinates.
    @param start_res: The number of the starting residue
    @param end_res: The number of the ending residue
    @return (coords, indeces): coords - A 3 x n matrix where n is the number of atoms in the matrix
        indeces - The indeces into this array, which indicate the position of the P atoms of each residue

    '''

    coords = []
    indeces = dict()
    count = 0

    for i in range(start_res, end_res+1):
        res = chain[i]
        indeces[res.id[1]] = count
        #print 'res:', res.resname.strip()
        for aname in a_names[res.resname.strip()]:
            coords += [res[aname].get_vector().get_array()]
            count += 1
            continue


    return (coords, indeces)

def get_atom_name_array(chain, start_res, end_res):
    '''
    Return an array of the coordinates of all the atoms in the chain, 
    arranged in the following order:

    P, O5', C5', C4', O4', C3', C1', C2', base_atoms, O3'

    @param chain: The chain from which to get the coordinates.
    @param start_res: The number of the starting residue
    @param end_res: The number of the ending residue
    @return (coords, indeces): coords - A 3 x n matrix where n is the number of atoms in the matrix
        indeces - The indeces into this array, which indicate the position of the P atoms of each residue

    '''

    coords = []
    indeces = dict()
    count = 0

    for i in range(start_res, end_res+1):
        res = chain[i]
        indeces[res.id[1]] = count
        #print 'res:', res.resname.strip()
        for aname in a_names[res.resname.strip()]:
            coords += ['%d%s' % (i, aname)]
            count += 1
            continue


    return (coords, indeces)

def set_atom_coord_array(chain, coords, start_res, end_res):
    '''
    Set the coordinates of the atoms in the chain to the ones in coords. 

    P, O5', C5', C4', O4', C3', C1', C2', base_atoms, O3'

    @param chain: The chain which will recieve coordinates
    @param coords: The coordinates to be entered
    @param start_res: The number of the starting residue
    @param end_res: The number of the ending residue
    @return (coords, indeces): coords - A 3 x n matrix where n is the number of atoms in the matrix
        indeces - The indeces into this array, which indicate the position of the P atoms of each residue

    '''
    count = 0

    for i in range(start_res, end_res+1):
        res = chain[i]
        #print 'res:', res.resname.strip()
        for aname in a_names[res.resname.strip()]:
            #print "coords[count]:", coords[count]
            #chain[i][aname].coord = bpdb.Vector(coords[count])
            chain[i][aname].coord = coords[count]
            #print "coords2:", chain[i][aname].coord
            count += 1
    return chain

def align_starts(chain_stems, chain_loop, handles, end=0):
    '''
    Align the sugar rings of one part of the stem structure to one part
    of the loop structure.

    @param chain_stems: The chain containing the stems
    @param chain_loop: The chain containing the sampled loop
    @param handles: The indexes into the stem and loop for the overlapping residues.
    '''
    if end == 0:
        v1 = get_alignment_vectors_rosetta(chain_stems, handles[0], handles[1])
        v2 = get_alignment_vectors_barnacle(chain_loop, handles[2], handles[3])
    else:
        v1 = get_measurement_vectors_rosetta(chain_stems, handles[0], handles[1])
        v2 = get_measurement_vectors_barnacle(chain_loop, handles[2], handles[3])

    #v1_m = get_measurement_vectors_rosetta(chain_stems, handles[0], handles[1])
    #v2_m = get_measurement_vectors_barnacle(chain_loop, handles[2], handles[3])

    v1_centroid = cuv.get_vector_centroid(v1)
    v2_centroid = cuv.get_vector_centroid(v2)

    v1_t = v1 - v1_centroid
    v2_t = v2 - v2_centroid

    sup = brmsd.optimal_superposition(v2_t, v1_t)

    #v2_mt = v2_m - v2_centroid
    #v2_rmt = np.dot(sup.transpose(), v2_mt.transpose()).transpose()

    #v1_rmt = v1_m - v1_centroid

    #rmsd = cuv.vector_set_rmsd(v1_rmt, v2_rmt)

    for atom in bpdb.Selection.unfold_entities(chain_loop, 'A'):
        atom.transform(np.eye(3,3), -v2_centroid)
        atom.transform(sup, v1_centroid)

def get_adjacent_interatom_distances(chain, start_res, end_res):
    adjacent_atoms = dict()
    adjacent_atoms['P'] = ['O5\'']
    adjacent_atoms['O5\''] = ['C5\'']
    adjacent_atoms['C5\''] = ['C4\'']
    adjacent_atoms['C4\''] = ['O4\'', 'C3\'']
    adjacent_atoms['O4\''] = ['C1\'']
    adjacent_atoms['C1\''] = ['C2\'']
    adjacent_atoms['C2\''] = ['C3\'']
    adjacent_atoms['C3\''] = ['O3\'']

    distances = []
    ress = list(chain.get_list())
    for i in range(start_res, end_res+1):
        res = chain[i]
        for key in adjacent_atoms.keys():
            for value in adjacent_atoms[key]:
                distances += [res[key] - res[value]]

    for i in range(start_res, end_res+1):
        #print "ress[i]['P'].coord:", ress[i]['P'].coord
        #print "ress[i]['O3\''].coord:", ress[i-1]['O3\''].coord
        distances += [ress[i]['P'] - ress[i-1]['O3\'']]

    return distances


def get_adjacent_interatom_names(chain, start_res, end_res):
    adjacent_atoms = dict()
    adjacent_atoms['P'] = ['O5\'']
    adjacent_atoms['O5\''] = ['C5\'']
    adjacent_atoms['C5\''] = ['C4\'']
    adjacent_atoms['C4\''] = ['O4\'', 'C3\'']
    adjacent_atoms['O4\''] = ['C1\'']
    adjacent_atoms['C1\''] = ['C2\'']
    adjacent_atoms['C2\''] = ['C3\'']
    adjacent_atoms['C3\''] = ['O3\'']

    distances = []
    ress = list(chain.get_list())
    for i in range(start_res, end_res+1):
        for key in adjacent_atoms.keys():
            for value in adjacent_atoms[key]:
                distances += [str(key) + "-" + str(value)]

    for i in range(start_res, end_res+1):
        #print "ress[i]['P'].coord:", ress[i]['P'].coord
        #print "ress[i]['O3\''].coord:", ress[i-1]['O3\''].coord
        distances += ['%dP-O3' % (i)]

    return distances


def close_fragment_loop(chain_stems, chain_loop, handles, iterations=1000):
    '''
    Align the chain_loop so that it stretches from the end of one stem to the 
    start of the other.

    @param chain_stems: The PDB coordinates of the chain with the stem.
    @param chain_loop: The PDB coordinates of the sampled loop.
    @param iterations: The number of iterations to use for the CCD loop closure.
    '''

    align_starts(chain_stems, chain_loop, handles)
    e = np.eye(3,3)

    (moving, indeces) = get_atom_coord_array(chain_loop, handles[2], handles[3])
    (moving_names, indeces_names) = get_atom_name_array(chain_loop, handles[2], handles[3])

    moving1 = np.array(moving)
    moving2 = np.array(moving)

    moving_orig = np.array(moving)

    moving = moving1
    tmoving = moving2

    fixed = np.array(get_measurement_vectors_rosetta(chain_stems, handles[0], handles[1]))

    start_res = handles[2]
    end_res = handles[3]

    start_index = indeces[handles[2]+1]
    end_index = indeces[handles[3]]

    points = []
    for i in range(handles[2]+1, handles[3]+1):
        si = indeces[i]
        points += [[si+1, si+2, si+3]]

    #points = [[start_index + 1, start_index + 2, start_index + 3], [end_index + 1, end_index + 2, end_index + 3]]
    #points = [[end_index + 1, end_index + 2, end_index + 3]]
    rot_mat = np.eye(3,3)
        
    #distances = get_adjacent_interatom_distances(chain_loop, handles[2], handles[3])
    #names = get_adjacent_interatom_names(chain_loop, handles[2], handles[3])

    for i in range(iterations):
        rmsd = cbc.calc_rmsd(moving[end_index+4:end_index+7], fixed)
        if rmsd < 0.1:
            break
        #print "rmsd:", rmsd
        #print "moving:", moving[:5]
        prev_l = points[0][0]

        for j in xrange(len(points)):
            for k in xrange(len(points[j])):
                l = points[j][k]
                #print "l:", l, "len(moving):", len(moving)
                #print "moving1:", moving[:6]
                point = moving[l-1]
                axis = moving[l] - moving[l-1]
                target = moving[end_index+4:end_index+7]
                #print moving_names[end_index+4:end_index+7]
                
                '''
                print "axis:", axis
                print "point:", point
                print "target:", target
                print "fixed:", fixed
                '''
                cbc.get_closer_rotation_matrix(axis, point, target, fixed, rot_mat)

                rem_moving = moving[l+1:]
                rem_tmoving = tmoving[l+1:]

                rem_moving -= moving[l]
                np.dot(rem_moving, rot_mat.transpose(), out = rem_tmoving)
                rem_tmoving += moving[l]

                #chain_loop = set_atom_coord_array(chain_loop, moving, handles[2], handles[3])
                #distances2 = get_adjacent_interatom_distances(chain_loop, handles[2], handles[3])

                '''
                print "i:", i, "l:", l
                print zip(moving_names, list(moving))
                d3 =  zip(names, list(np.array(distances) - np.array(distances2)))
                for z in range(len(d3)):
                    #print "z1:", d3[z][1]
                    if abs(d3[z][1]) > 0.1:
                        print d3[z]

                assert(np.allclose(distances, distances2))
                '''
                tt_moving = moving
                moving = tmoving

                moving[prev_l+1:l+1] = tt_moving[prev_l+1:l+1]
                tmoving = tt_moving

                prev_l = l

    #print "moving:", moving
    
    assert(not np.allclose(moving_orig, moving))

    chain_loop = set_atom_coord_array(chain_loop, moving, handles[2], handles[3])
    (moving_new, indeces) = get_atom_coord_array(chain_loop, handles[2], handles[3])

    assert(not np.allclose(moving_orig, moving_new))

    return (rmsd, chain_loop)

def quality_measurement(chain_stems, chain_barnacle, handles):
    '''
    See how well we can fit the ends of the loop to the ends of the
    stems.

    First the suger ring of one of the residues in the loop (barnacle) chain
    is superimposed onto the suger ring of the end of the stem residue.

    The torsion angles of the first and last residues of the loop region
    are then modified to try and close the loop between the one end
    of the first stem and the other end of the other stem (which can also
    be part of the same stem in the case of a loop region).

    @param chain_stems: The chain containing the reconstructed stems.
    @param chain_barnacle: The chain which contains the sampled loop region.
    @param handles: The residues of the ends of the stems as well as the 
        indeces into the loop chain.
    '''
    align_starts(chain_stems, chain_barnacle, handles)

    v1_m = get_measurement_vectors_rosetta(chain_stems, handles[0], handles[1])
    v2_m = get_measurement_vectors_barnacle(chain_barnacle, handles[2], handles[3])

    rmsd2 = get_measurement_rmsd(chain_stems, chain_barnacle, handles)

    for i in range(50):
        perturb_c3_o3(chain_barnacle, handles[2], handles[3], v1_m)

        for j in range(handles[2] + 1, handles[2] + 2):
            perturb_p_o5(chain_barnacle, j, handles[3], v1_m)
            perturb_o5_c5(chain_barnacle, j, handles[3], v1_m)
            perturb_c5_c4(chain_barnacle,  j, handles[3], v1_m)

        perturb_p_o5(chain_barnacle, handles[3], handles[3], v1_m)
        perturb_o5_c5(chain_barnacle, handles[3], handles[3], v1_m)
        perturb_c5_c4(chain_barnacle, handles[3], handles[3], v1_m)


    rmsd9 = get_measurement_rmsd(chain_stems, chain_barnacle, handles)
    return rmsd9, chain_barnacle

def trim_chain(chain, start_res, end_res):
    '''
    Remove all residues that are not between start_res and end_res.
    '''
    to_detach = []
    for res in chain:
        if res.id[1] < start_res+1 or end_res-1 < res.id[1]:
            to_detach += [res]

    for res in to_detach:
        chain.detach_child(res.id)

def reconstruct_loop(chain, sm, ld, side=0):
    '''
    Reconstruct a particular loop.

    The chain should already have the stems reconstructed.

    @param chain: A Bio.PDB.Chain structure.
    @param sm: A SpatialModel structure
    @param ld: The name of the loop
    '''

    bg = sm.bg
    seq = bg.get_flanking_sequence(ld, side)
    (a,b,i1,i2) = bg.get_flanking_handles(ld, side)

    model = barn.Barnacle(seq)
    model.sample()
    s = model.structure

    prev_r = 1000.
    min_r = 1000.
    min_contacts = (1000, 100.)
    iterations = 20
    best_loop_chain = None
    sys.stdout.write("reconstructing %s:" % (ld))

    for i in range(iterations):
        model.sample()
        chain_loop = list(model.structure.get_chains())[0]

        loop_atoms = bpdb.Selection.unfold_entities(chain_loop, 'A') 
        
        ns = bpdb.NeighborSearch(loop_atoms)
        contacts1 = len(ns.search_all(0.8))
        
        if a == 0:
            align_starts(chain, chain_loop, (a,b,i1,i2), end=1)
            loop_chain = chain_loop
            r = 0.000
        elif b == bg.length:
            align_starts(chain, chain_loop, (a,b,i1,i2), end=0)
            loop_chain = chain_loop
            r = 0.000
        else:
            r, loop_chain = close_fragment_loop(chain, chain_loop, (a,b,i1,i2), iterations=2000)

        '''
        sample_len = ss.poisson.rvs(2)

        if sample_len == 0:
            sample_len = 1

        j = rand.randint(0, len(seq)-sample_len)
        m1 = j
        m2 = m1 + sample_len

        model.sample(start=m1, end=m2)

        # attempt to close the distance between the two stems with this loop and return
        # the minimum distance achieved
        r, loop_chain = close_fragment_loop(chain, list(model.structure.get_chains())[0], (a,b,i1,i2))

        multiplier = .001 ** (1 / float(iterations))
        temperature = 1. * (multiplier) ** i

        if r > prev_r:
            factor = -(r - prev_r) / ( temperature)
            
            p = np.exp (factor)
            
            if rand.random() > p:
                model.undo()
                continue

        prev_r = r
        '''
        orig_loop_chain = copy.deepcopy(loop_chain)

        trim_chain(loop_chain, i1, i2)
        loop_atoms = bpdb.Selection.unfold_entities(loop_chain, 'A')

        if a != 0:
            loop_atoms += bpdb.Selection.unfold_entities(chain[a], 'A')
        if b != bg.length:
            loop_atoms += bpdb.Selection.unfold_entities(chain[b], 'A')

        if a != 0:
            loop_atoms += bpdb.Selection.unfold_entities(chain[a-1], 'A')
        if b != bg.length:
            loop_atoms += bpdb.Selection.unfold_entities(chain[b+1], 'A')

        ns = bpdb.NeighborSearch(loop_atoms)

        contacts2 = len(ns.search_all(1.0))
        sys.stdout.write('.')
        sys.stdout.flush()

        #print "r:", r, "contacts1:", contacts1, "contacts2:",  contacts2

        if (contacts2, r) < min_contacts:
            best_loop_chain = copy.deepcopy(orig_loop_chain)
            min_contacts = (contacts2, r)
            #print "min_contacts:", min_contacts

        '''
        if contacts2 == 0:
            break
        '''

        trim_chain(loop_chain, i1, i2)
        output_chain(loop_chain, os.path.join(conf.Configuration.test_output_dir, 's%d.pdb' % (i+2)))
        
    #output_chain(chain, os.path.join(conf.Configuration.test_output_dir, 's1.pdb'))
    trim_chain(best_loop_chain, i1, i2+1)
    output_chain(best_loop_chain, os.path.join(conf.Configuration.test_output_dir, 'sb.pdb'))
    add_loop_chain(chain, best_loop_chain, (a,b,i1,i2), bg.length)
    sys.stdout.write('\n')
    sys.stdout.flush()

def reconstruct_loops(chain, sm):
    '''
    Reconstruct the loops of a chain.

    All of the stems should already be reconstructed in chain.

    @param chain: A Bio.PDB.Chain chain.
    @param sm: The SpatialModel from which to reconstruct the loops.
    '''
    for d in sm.bg.defines.keys():
        if d[0] != 's':
            if sm.bg.weights[d] == 2:
                reconstruct_loop(chain, sm, d, 0)
                reconstruct_loop(chain, sm, d, 1)
            else:
                reconstruct_loop(chain, sm, d)
            

def reconstruct(sm):
    '''
    Re-construct a full-atom model from a coarse-grain model.

    @param bg: The BulgeGraph
    @return: A Bio.PDB.Chain chain
    '''
    chain = reconstruct_stems(sm)
    reconstruct_loops(chain, sm)
    return chain
