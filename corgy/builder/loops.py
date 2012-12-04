import random as rand
import numpy as np
import copy

import Bio.PDB as bpdb

import corgy.aux.ccd.cytvec as cv
import corgy.aux.CPDB.src.examples.BarnacleCPDB as barn

import corgy.builder.ccd as cbc
import corgy.builder.rmsd as brmsd

import corgy.utilities.debug as cud
import corgy.utilities.pdb as cup
import corgy.utilities.vector as cuv

import sys

class BarnacleEnergy:
    def __init__(self):
        self.prev_r = 1000.
        self.min_r = 1000.
        self.min_contacts = 1000.
    
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

def get_alignment_vectors(ress, r1, r2):
    return( ress[r1]['C4*'].get_vector().get_array(),
            ress[r1]['C3*'].get_vector().get_array(),
            ress[r1]['O3*'].get_vector().get_array())

def get_measurement_vectors(ress, r1, r2):
    #cud.pv('bpdb.Selection.unfold_entities(ress[r2], \'A\')')
    return( ress[r2]['C2*'].get_vector().get_array(), 
            ress[r2]['C3*'].get_vector().get_array(),
            ress[r2]['O3*'].get_vector().get_array())

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
            try:
                coords += [res[aname].get_vector().get_array()]
            except KeyError:
                #alternate name for the OP1 atoms
                if aname == 'OP1':
                    coords += [res['O1P'].get_vector().get_array()]
                elif aname == 'OP2':
                    coords += [res['O2P'].get_vector().get_array()]
                else:
                    raise
                    
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
            try:
                chain[i][aname].coord = coords[count]
            except KeyError:
                if aname == 'OP1':
                    chain[i]['O1P'].coord = coords[count]
                elif aname == 'OP2':
                    chain[i]['O2P'].coord = coords[count]
                else:
                    raise
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
    #(moving_names, indeces_names) = get_atom_name_array(chain_loop, handles[2], handles[3])


    fixed = np.array(get_measurement_vectors(chain_stems, handles[0], handles[1]))

    start_res = handles[2]
    end_res = handles[3]

    start_index = indeces[handles[2]+1]
    end_index = indeces[handles[3]+1]

    points = []
    #points += [indeces[handles[2]+1]]

    #points += indeces[handles[2]+1] #O3* -> P bond
    #for i in range(handles[2]+1, handles[3]+1):
    for i in [handles[2]+1, handles[3]]:
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

def align_and_close_loop(seq_len, chain, chain_loop, (a,b,i1,i2)):
    '''
    Align chain_loop to the scaffold present in chain.

    This means that nt i1 in chain_loop will be aligned
    to nt a, and nt i2 will attempt to be aligned to nt b.

    @param seq_len: The length of the entire sequence (including stems)
    @param chain: The scaffold containing the stems
    @param chain_loop: The sampled loop region

    @para (a,b,i1,i2): The handles indicating which nucleotides
                       will be aligned

    @return: (r, loop_chain), where r is the rmsd of the closed loop
            and loop_chain is the closed loop
    '''
    loop_atoms = bpdb.Selection.unfold_entities(chain_loop, 'A') 
    
    ns = bpdb.NeighborSearch(loop_atoms)
    contacts1 = len(ns.search_all(0.8))
    
    if a == 0:
        align_starts(chain, chain_loop, (a,b,i1,i2), end=1)
        loop_chain = chain_loop
        r = 0.000
    elif b == seq_len:
        align_starts(chain, chain_loop, (a,b,i1,i2), end=0)
        loop_chain = chain_loop
        r = 0.000
    else:
        r, loop_chain = close_fragment_loop(chain, chain_loop, (a,b,i1,i2), iterations=10000)

    return (r, loop_chain)

def build_loop(stem_chain, loop_seq, (a,b,i1,i2), seq_len, iterations, consider_contacts=False):
    '''
    Build a loop.

    The alignment vectors of the nucleotides a and b in stem_chain
    should align to the alignment vectors of nucleotides i1 and i2
    of the sequence sampled by loop_seq.

    @param stem_chain: A chain containing the assembled stems.
    @param loop_seq: The sequence of the loop region including the
                     adjacent stem segments
    @param (a,b,i1,i2): The numbers of the nucleotides defining
                     where the loop starts in the whole sequence (a and b)
                     and within the loop sequence (i1 and i2)
    @param seq_len: The length of the entire sequence (including the stems)
    @param iterations: The number of MCMC iterations to run
    @param consider_contacts: The into account the contact between the loop and other portions
                              of the structure

    @return: A Bio.PDB.Chain structure containing the best sampled loop.
    '''
    model = barn.BarnacleCPDB(loop_seq, 2.)

    best_loop_chain = None
    prev_r = 1000.
    min_r = 1000.
    min_contacts = (10000000., 1000.)
    cud.pv('cup.num_noncovalent_clashes(stem_chain)')

    for i in range(iterations):
        start = rand.randint(0, len(loop_seq)-1)
        end = rand.randint(start+1, len(loop_seq))

        #print "start:", start, "end:", end
        #model.sample(start=start, end=end)

        model.sample()
        chain_loop = list(model.structure.get_chains())[0]

        (r, loop_chain) = align_and_close_loop(seq_len, stem_chain, chain_loop, (a, b, i1, i2))

        #orig_loop_chain = copy.deepcopy(loop_chain)
        orig_loop_chain = copy.deepcopy(loop_chain)

        all_chain = copy.deepcopy(stem_chain)
        cup.trim_chain(loop_chain, i1, i2+1)
        add_loop_chain(all_chain, loop_chain, (a,b,i1,i2), seq_len)
        contacts2 = cup.num_noncovalent_clashes(all_chain)

        sys.stderr.write('.')
        sys.stderr.flush()

        #print "r:", r, "contacts1:", contacts1, "contacts2:",  contacts2
        if consider_contacts: 
            if (contacts2, r) < min_contacts:
            #if (0, r) < min_contacts:
                best_loop_chain = copy.deepcopy(orig_loop_chain)
                min_contacts = (contacts2, r)
                #min_contacts = (0, r)
                print "min_contacts:", min_contacts
        else:
            if (0, r) < min_contacts:
                best_loop_chain = copy.deepcopy(orig_loop_chain)
                min_contacts = (0, r)
                #print "r:", r
        '''
        if contacts2 == 0:
            break
        '''

        #trim_chain(loop_chain, i1, i2)
    sys.stderr.write(str(min_contacts))
    return best_loop_chain
