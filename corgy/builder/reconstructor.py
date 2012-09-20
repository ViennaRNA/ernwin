import corgy.builder.models as models
import corgy.builder.rmsd as brmsd

import corgy.graph.graph_pdb as gpdb
import corgy.utilities.vector as cuv

import corgy.builder.ccd as cbc
import corgy.aux.ccd.cytvec as cv

import corgy.aux.Barnacle as barn

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

def get_alignment_vectors(ress, r1, r2):
    return( ress[r1]['C4*'].get_vector().get_array(),
            ress[r1]['C3*'].get_vector().get_array(),
            ress[r1]['O3*'].get_vector().get_array())

def get_measurement_vectors1(ress, r1, r2):
    return( ress[r2]['C4*'].get_vector().get_array(), 
            ress[r2]['C3*'].get_vector().get_array(),
            ress[r2]['O3*'].get_vector().get_array())

def get_measurement_vectors2(ress, r1, r2):
    return( ress[r2]['O4*'].get_vector().get_array(), 
            ress[r2]['C1*'].get_vector().get_array(),
            ress[r2]['C2*'].get_vector().get_array())

def get_measurement_vectors(ress, r1, r2):
    return( ress[r2]['C2*'].get_vector().get_array(), 
            ress[r2]['C3*'].get_vector().get_array(),
            ress[r2]['O3*'].get_vector().get_array())

def pdb_rmsd(c1, c2):
    '''
    Calculate the all-atom rmsd between two RNA chains.

    @param c1: A Bio.PDB.Chain
    @param c2: Another Bio.PDB.Chain
    @return: The rmsd between the locations of all the atoms in the chains.
    '''

    a_5_names = ['P', 'O5*', 'C5*', 'C4*', 'O4*', 'O2*']
    a_3_names = ['C1*', 'C2*', 'C3*', 'O3*']

    a_names = dict()
    a_names['U'] = a_5_names + ['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6'] + a_3_names
    a_names['C'] = a_5_names + ['N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6'] + a_3_names

    a_names['A'] = a_5_names + ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N6', 'N7', 'C8', 'N9'] + a_3_names
    a_names['G'] = a_5_names + ['N1', 'C2', 'N2', 'N3', 'C4', 'C5', 'C6', 'O6', 'N7', 'C8', 'N9'] + a_3_names

    all_atoms1 = []
    all_atoms2 = []

    if len(c1.get_list()) != len(c2.get_list()):
        print >>sys.stderr, "Chains of different length"
        raise Exception("Chains of different length.")

    for i in range(1, len(list(c1.get_list()))+1):
        anames = a_5_names + a_names[c1[i].resname.strip()] + a_3_names
        #anames = a_5_names + a_3_names

        try:
            atoms1 = [c1[i][a] for a in anames]
            atoms2 = [c2[i][a] for a in anames]
        except KeyError:
            print >>sys.stderr, "Residue number %d is missing an atom, continuing with the rest." % (i)
            continue

        if len(atoms1) != len(atoms2):
            print >>sys.stderr, "Number of atoms differs in the two chains."
            raise Exception("Missing atoms.")

        all_atoms1 += atoms1
        all_atoms2 += atoms2

    print "rmsd len:", len(all_atoms1), len(all_atoms2)
    sup = bpdb.Superimposer()
    sup.set_atoms(all_atoms1, all_atoms2)

    return sup.rms

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
        atom.transform(cuv.identity_matrix, translation)

def align_chain_to_stem(chain, define, stem2):
    stem1 = define_to_stem_model(chain, define)

    (r, u, v, t) = gpdb.get_stem_orientation_parameters(stem1.vec(), stem1.twists[0], stem2.vec(), stem2.twists[0])
    rot_mat = get_stem_rotation_matrix(stem1, (math.pi-u, -v, -t))
    rotate_chain(chain, np.linalg.inv(rot_mat), stem1.mids[0])
    translate_chain(chain, stem2.mids[0] - stem1.mids[0])

def reconstruct_stems(sm, stem_library=dict()):
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
        #print "stem_name:", stem_name, "stem_def:", stem_def, "filename:", filename
        pdb_file = os.path.join(conf.Configuration.stem_fragment_dir, filename)

        #print len(stem_library.keys())
        if filename in stem_library.keys():
            chain = stem_library[filename]
        else:
            chain = list(bpdb.PDBParser().get_structure('temp', pdb_file).get_chains())[0]
            stem_library[filename] = chain

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
        loop_chain[handles[3]].id = (' ', handles[1], ' ')
        add_residue_to_rosetta_chain(chain, loop_chain[handles[3]])
    else:
        loop_chain[handles[3]].id = (' ', handles[1], ' ')
        add_residue_to_rosetta_chain(chain, loop_chain[handles[3]])

    # We won't replace the last residue
    counter = 1
    for i in range(handles[2]+1, handles[3]):
        loop_chain[i].id = (' ', handles[0] + counter, ' ')
        add_residue_to_rosetta_chain(chain, loop_chain[i])
        counter += 1


a_5_names = ['P', 'OP1', 'OP2', 'P', 'O5*', 'C5*', 'C4*', 'O4*', 'O2*']
a_3_names = ['C1*', 'C2*', 'C3*', 'O3*']

side_chain_atoms = dict()

side_chain_atoms['U'] = ['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6']
side_chain_atoms['C'] = ['N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6']

side_chain_atoms['A'] = ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N6', 'N7', 'C8', 'N9']
side_chain_atoms['G'] = ['N1', 'C2', 'N2', 'N3', 'C4', 'C5', 'C6', 'O6', 'N7', 'C8', 'N9']

a_names = dict()
a_names['U'] = a_5_names + side_chain_atoms['U'] + a_3_names
a_names['C'] = a_5_names + side_chain_atoms['C'] + a_3_names

a_names['A'] = a_5_names + side_chain_atoms['A'] + a_3_names
a_names['G'] = a_5_names + side_chain_atoms['G'] + a_3_names

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

    for i in range(start_res, end_res+2):
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

    for i in range(start_res, end_res+2):
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

    for i in range(start_res, end_res+2):
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
        v1 = get_alignment_vectors(chain_stems, handles[0], handles[1])
        v2 = get_alignment_vectors(chain_loop, handles[2], handles[3])
    else:
        v1 = get_measurement_vectors(chain_stems, handles[0], handles[1])
        v2 = get_measurement_vectors(chain_loop, handles[2], handles[3])

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
    adjacent_atoms['P'] = ['O5*', 'OP1', 'OP2']
    adjacent_atoms['O5*'] = ['C5*']
    adjacent_atoms['C5*'] = ['C4*']
    adjacent_atoms['C4*'] = ['O4*', 'C3*']
    adjacent_atoms['O4*'] = ['C1*']
    adjacent_atoms['C1*'] = ['C2*']
    adjacent_atoms['C2*'] = ['C3*', 'O2*']
    adjacent_atoms['C3*'] = ['O3*']

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
        distances += [ress[i]['P'] - ress[i-1]['O3*']]

    return distances


def get_adjacent_interatom_names(chain, start_res, end_res):
    adjacent_atoms = dict()
    adjacent_atoms['P'] = ['O5*', 'OP1', 'OP2']
    adjacent_atoms['O5*'] = ['C5*']
    adjacent_atoms['C5*'] = ['C4*']
    adjacent_atoms['C4*'] = ['O4*', 'C3*']
    adjacent_atoms['O4*'] = ['C1*']
    adjacent_atoms['C1*'] = ['C2*']
    adjacent_atoms['C2*'] = ['C3*', 'O2*']
    adjacent_atoms['C3*'] = ['O3*']

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


def close_fragment_loop(chain_stems, chain_loop, handles, iterations=5000):
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


    fixed = np.array(get_measurement_vectors(chain_stems, handles[0], handles[1]))

    start_res = handles[2]
    end_res = handles[3]

    start_index = indeces[handles[2]+1]
    end_index = indeces[handles[3]+1]

    points = []
    #points += [indeces[handles[2]+1]]

    #points += indeces[handles[2]+1] #O3* -> P bond
    for i in range(handles[2]+1, handles[3]+1):
        si = indeces[i]

        # 
        points += [si]
        points += [si+4, si+5, si+6]

    rot_mat = np.eye(3,3)
        
    '''
    distances = get_adjacent_interatom_distances(chain_loop, handles[2], handles[3])
    names = get_adjacent_interatom_names(chain_loop, handles[2], handles[3])
    '''


    moving = np.array(moving)
    points = np.array(points)

    cv.ccd_cython(moving, fixed, points, end_index-3, iterations) 
    rmsd = cbc.calc_rmsd(moving[end_index-3:end_index], fixed)

    chain_loop = set_atom_coord_array(chain_loop, moving, handles[2], handles[3])

    '''
    assert(not np.allclose(moving_orig, moving))

    (moving_new, indeces) = get_atom_coord_array(chain_loop, handles[2], handles[3])

    assert(not np.allclose(moving_orig, moving_new))

    distances2 = get_adjacent_interatom_distances(chain_loop, handles[2], handles[3])
    print [ (d,n) for (d,n) in zip(np.array(distances2) - np.array(distances), names) if (d > 0.0001 or math.isnan(d)) ]

    #print np.array(distances2) - np.array(distances)

    assert(np.allclose(distances, distances2))
    '''


    return (rmsd, chain_loop)

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

def print_alignment_pymol_file(handles):
    output_str = """
select bb, /s2///%d/O4* | /s2///%d/C1* | /s2///%d/C1*
show sticks, bb
color red, bb

select bb, /s2///%d/O4* | /s2///%d/C1*| /s2///%d/C2*
show sticks, bb
color red, bb

select bb, s1///%d/O4* | s1///%d/C1* | s1///%d/C2*
show sticks, bb
color green, bb

select bb, s1///%d/O4* | s1///%d/C1* | s1///%d/C2*
show sticks, bb
color green, bb

show cartoon, all
""" % (handles[2], handles[2], handles[2],
        handles[3], handles[3], handles[3],
        handles[0], handles[0], handles[0],
        handles[1], handles[1], handles[1])
    output_file = os.path.join(conf.Configuration.test_output_dir, "align.pml")
    f = open(output_file, 'w')
    f.write(output_str)
    f.flush()
    f.close()


def reconstruct_loop(chain, sm, ld, side=0, samples=40):
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

    #print "ld:", ld, "(a,b,i1,i2)", a,b,i1,i2
    #print "seq:", seq

    model = barn.Barnacle(seq)
    model.sample()
    s = model.structure

    prev_r = 1000.
    min_r = 1000.
    min_contacts = (1000, 100.)
    iterations = samples
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
            r, loop_chain = close_fragment_loop(chain, chain_loop, (a,b,i1,i2), iterations=500)

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
        #orig_loop_chain = copy.deepcopy(loop_chain)
        orig_loop_chain = loop_chain

        trim_chain(loop_chain, i1, i2+1)
        loop_atoms = bpdb.Selection.unfold_entities(loop_chain, 'A')
        #loop_atoms += bpdb.Selection.unfold_entities(chain, 'A')

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

        #trim_chain(loop_chain, i1, i2)
        
    sys.stdout.write(str(min_contacts))
    output_chain(chain, os.path.join(conf.Configuration.test_output_dir, 's1.pdb'))
    output_chain(best_loop_chain, os.path.join(conf.Configuration.test_output_dir, 's2.pdb'))
    print_alignment_pymol_file((a,b,i1,i2))

    trim_chain(best_loop_chain, i1, i2+1)
    add_loop_chain(chain, best_loop_chain, (a,b,i1,i2), bg.length)
    sys.stdout.write('\n')
    sys.stdout.flush()

def reconstruct_loops(chain, sm, samples=40):
    '''
    Reconstruct the loops of a chain.

    All of the stems should already be reconstructed in chain.

    @param chain: A Bio.PDB.Chain chain.
    @param sm: The SpatialModel from which to reconstruct the loops.
    '''
    for d in sm.bg.defines.keys():
        if d[0] != 's':
            if sm.bg.weights[d] == 2:
                reconstruct_loop(chain, sm, d, 0, samples=samples)
                reconstruct_loop(chain, sm, d, 1, samples=samples)
            else:
                reconstruct_loop(chain, sm, d, 0, samples=samples)
            

def reconstruct(sm):
    '''
    Re-construct a full-atom model from a coarse-grain model.

    @param bg: The BulgeGraph
    @return: A Bio.PDB.Chain chain
    '''
    chain = reconstruct_stems(sm)
    reconstruct_loops(chain, sm)
    return chain

def replace_base(res_dir, res_ref):
    '''
    Orient res_ref so that it points in the same direction
    as res_dir.

    @param res_dir: The residue indicating the direction
    @param res_ref: The reference residue to be rotated
    @return res: A residue with the atoms of res_ref pointing in the direction of res_dir
    '''
    #av = { 'U': ['N1', 'C6', 'C2'], 'C': ['N1', 'C6', 'C2'], 'A': ['N9', 'C4', 'C8'], 'G': ['N9', 'C4', 'C8'] }
    av = { 'U': ['N1', 'C1*', 'C2*'], 'C': ['N1', 'C1*', 'C2*'], 'A': ['N9', 'C1*', 'C2*'], 'G': ['N9', 'C1*', 'C2*'] }

    dv = av[res_dir.resname.strip()]
    rv = av[res_ref.resname.strip()]

    dir_points = np.array([res_dir[v].get_vector().get_array() for v in dv])
    ref_points = np.array([res_ref[v].get_vector().get_array() for v in rv])

    dir_centroid = cuv.get_vector_centroid(dir_points)
    ref_centroid = cuv.get_vector_centroid(ref_points)

    #sup = brmsd.optimal_superposition(dir_points - dir_centroid, ref_points - ref_centroid)
    sup = brmsd.optimal_superposition(ref_points - ref_centroid, dir_points - dir_centroid)
    new_res = copy.deepcopy(res_ref)

    for atom in new_res:
        atom.transform(np.eye(3,3), -ref_centroid)
        atom.transform(sup, dir_centroid)

    #print "dir_points:", dir_points
    #print "ref_points:", ref_points
    return new_res

def replace_bases(chain, seq):
    '''
    Go through the chain and replace the bases with the ones specified in the
    sequence.

    This is necessary since the stems are sampled by their length rather than
    sequence. Therefore some stem fragments may contain a sequence that is
    different from the one that is required.

    This method will change the identity of those bases, as well as align the atoms
    according to their N1->C2, N1->C6 or N9->C4, N9->C8. vector pairs.

    param @chain: A Bio.PDB.Chain with some reconstructed residues
    param @seq: The sequence of the structure represented by chain
    '''
    s1 = bpdb.PDBParser().get_structure('t', conf.Configuration.template_residue_fn)
    tchain = list(s1.get_chains())[0]


    tindeces = { 'A': 1, 'C': 2, 'G': 3, 'U': 4}

    print "len(seq):", len(seq)
    ress = chain.get_list()

    for i in range(len(ress)):
        num = ress[i].id[1]
        name = ress[i].resname.strip()

        ref_res = tchain[tindeces[seq[num-1]]]
        new_res = replace_base(ress[i], ref_res)

        sca = side_chain_atoms[ress[i].resname.strip()]
        for aname in sca:
            ress[i].detach_child(aname)
        
        sca = side_chain_atoms[new_res.resname.strip()]
        for aname in sca:
            ress[i].add(new_res[aname])

        ress[i].resname = new_res.resname
        '''
        ress[i].resname = new_res.resname
        ress[i].child_list = new_res.child_list
        ress[i].child_dict = new_res.child_dict
        '''

        
