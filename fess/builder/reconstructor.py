import multiprocessing as mp
import warnings

import itertools as it
import fess.builder.models as models
import os.path as op

import forgi.threedee.utilities.pdb as ftup
import forgi.threedee.utilities.average_atom_positions as ftua
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.vector as cuv
import forgi.utilities.debug as fud

import fess.builder.ccd as cbc

import borgy.builder.rmsd as brmsd

#import borgy.aux.Barnacle as barn
#import fess.aux.CPDB.BarnacleCPDB as barn

import Bio.PDB as bpdb
import Bio.PDB.Chain as bpdbc
import Bio.PDB.Model as bpdbm
import Bio.PDB.Structure as bpdbs

import scipy.stats as ss

from scipy.stats import norm, poisson

import os, math, sys
import borgy.builder.config as conf
import copy, time
import random as rand

import numpy as np

def get_measurement_vectors1(ress, r1, r2):
    return( ress[r2]["C4'"].get_vector().get_array(), 
            ress[r2]["C3'"].get_vector().get_array(),
            ress[r2]["O3'"].get_vector().get_array())

def get_measurement_vectors2(ress, r1, r2):
    return( ress[r2]["O4'"].get_vector().get_array(), 
            ress[r2]["C1'"].get_vector().get_array(),
            ress[r2]["C2'"].get_vector().get_array())

def rotate_stem(stem, (u, v, t)):
    '''
    Rotate a particular stem.
    '''
    stem2 = copy.deepcopy(stem)
    rot_mat4 = models.get_stem_rotation_matrix(stem, (u,v,t))
    stem2.rotate(rot_mat4, offset=stem.mids[0])

    return stem2

def reconstruct_stems(sm, stem_library=dict()):
    '''
    Reconstruct the stems around a Spatial Model.

    @param sm: Spatial Model
    '''
    #sm.traverse_and_build()
    new_chain = bpdbc.Chain(' ')

    for stem_name in sm.stem_defs.keys():
        models.reconstruct_stem(sm, stem_name, new_chain, stem_library)

    return new_chain

def output_chain(chain, filename, fr=None, to=None):
    '''
    Dump a chain to an output file. Remove the hydrogen atoms.

    @param chain: The Bio.PDB.Chain to dump.
    @param filename: The place to dump it.
    '''
    class HSelect(bpdb.Select):
        def accept_atom(self, atom):
            if atom.name.find('H') >= 0:
                return False
            else:
                return True

    m = bpdbm.Model(' ')
    s = bpdbs.Structure(' ')

    m.add(chain)
    s.add(m)

    io = bpdb.PDBIO()
    io.set_structure(s)
    io.save(filename, HSelect())

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

def print_alignment_pymol_file(handles):
    output_str = """
select bb, /s2///%d/O4' | /s2///%d/C1' | /s2///%d/C1'
show sticks, bb
color red, bb

select bb, /s2///%d/O4' | /s2///%d/C1'| /s2///%d/C2'
show sticks, bb
color red, bb

select bb, s1///%d/O4' | s1///%d/C1' | s1///%d/C2'
show sticks, bb
color green, bb

select bb, s1///%d/O4' | s1///%d/C1' | s1///%d/C2'
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

def get_flanking_stem_vres_distance(bg, ld):
    '''
    Get the distance between the two virtual residues adjacent
    to this bulge region.

    @param bg: The BulgeGraph data structure
    @param ld: The name of the linking bulge
    '''

    if len(bg.edges[ld]) == 2:
        connecting_stems = list(bg.edges[ld])

        (s1b, s1e) = bg.get_sides(connecting_stems[0], ld)
        (s2b, s2e) = bg.get_sides(connecting_stems[1], ld)

        if s1b == 1:
            (vr1_p, vr1_v, vr1_v_l, vr1_v_r) = ftug.virtual_res_3d_pos(bg, connecting_stems[0], bg.stem_length(connecting_stems[0]) - 1)
        else:
            (vr1_p, vr1_v, vr1_v_l, vr1_v_r) = ftug.virtual_res_3d_pos(bg, connecting_stems[0], 0)

        if s2b == 1:
            (vr2_p, vr2_v, vr2_v_l, vr2_v_r) = ftug.virtual_res_3d_pos(bg, connecting_stems[1], bg.stem_length(connecting_stems[1]) - 1)
        else:
            (vr2_p, vr2_v, vr2_v_l, vr2_v_r) = ftug.virtual_res_3d_pos(bg, connecting_stems[1], 0)

        dist2 = cuv.vec_distance((vr1_p + 7 * vr1_v), (vr2_p + 7. * vr2_v))
    else:
        dist2 = 0.

    return dist2

a_5_names = ["P", "OP1", "OP2", "P", "O5'", "C5'", "C4'", "O4'", "O2'"]
a_3_names = ["C1'", "C2'", "C3'", "O3'"]

backbone_atoms = ['P', "O5'", "C5'", "C4'", "C3'", "O3'"]

side_chain_atoms = dict()

side_chain_atoms['U'] = ['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6']
side_chain_atoms['C'] = ['N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6']

side_chain_atoms['A'] = ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N6', 'N7', 'C8', 'N9']
side_chain_atoms['G'] = ['N1', 'C2', 'N2', 'N3', 'C4', 'C5', 'C6', 'O6', 'N7', 'C8', 'N9']
side_chain_atoms['rU'] = ['N1', 'C2', 'O2', 'N3', 'C4', 'O4', 'C5', 'C6']
side_chain_atoms['rC'] = ['N1', 'C2', 'O2', 'N3', 'C4', 'N4', 'C5', 'C6']

side_chain_atoms['rA'] = ['N1', 'C2', 'N3', 'C4', 'C5', 'C6', 'N6', 'N7', 'C8', 'N9']
side_chain_atoms['rG'] = ['N1', 'C2', 'N2', 'N3', 'C4', 'C5', 'C6', 'O6', 'N7', 'C8', 'N9']

a_names = dict()
a_names['U'] = a_5_names + side_chain_atoms['U'] + a_3_names
a_names['C'] = a_5_names + side_chain_atoms['C'] + a_3_names

a_names['A'] = a_5_names + side_chain_atoms['A'] + a_3_names
a_names['G'] = a_5_names + side_chain_atoms['G'] + a_3_names

a_names['rU'] = a_5_names + side_chain_atoms['U'] + a_3_names
a_names['rC'] = a_5_names + side_chain_atoms['C'] + a_3_names

a_names['rA'] = a_5_names + side_chain_atoms['A'] + a_3_names
a_names['rG'] = a_5_names + side_chain_atoms['G'] + a_3_names

def get_alignment_vectors(ress, r1, r2):
    '''
    return( ress[r1]["C4'"].get_vector().get_array(),
            ress[r1]["C3'"].get_vector().get_array(),
            ress[r1]["O3'"].get_vector().get_array())
    '''
    vec = []
    for aname in backbone_atoms:
        try:
            r1_v = [ress[r1][aname].get_vector().get_array()]
            r2_v = [ress[r2][aname].get_vector().get_array()]

            vec += r1_v
            vec += r2_v
        except KeyError as ke:
            raise ke

    return vec

def get_measurement_vectors(ress, r1, r2):
    return( ress[r2]["C2'"].get_vector().get_array(), 
            ress[r2]["C3'"].get_vector().get_array(),
            ress[r2]["O3'"].get_vector().get_array())

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
        for aname in a_names[res.resname.strip()]:
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

    v1 = []
    v2 = []

    for handle in handles:
        if end == 0:
            v1 += get_alignment_vectors(chain_stems, handle[0], handle[1])
            v2 += get_alignment_vectors(chain_loop, handle[2], handle[3])
        elif end == 2:
            v1 = (chain_stems[handle[0]]["C4'"].get_vector().get_array(),
                  chain_stems[handle[0]]["C3'"].get_vector().get_array(),
                  chain_stems[handle[0]]["O3'"].get_vector().get_array())
            v2 = (chain_loop[handle[2]]["C4'"].get_vector().get_array(),
                  chain_loop[handle[2]]["C3'"].get_vector().get_array(),
                  chain_loop[handle[2]]["O3'"].get_vector().get_array())
        else:
            v1 += get_measurement_vectors(chain_stems, handle[0], handle[1])
            v2 += get_measurement_vectors(chain_loop, handle[2], handle[3])

    #v1_m = get_measurement_vectors_rosetta(chain_stems, handle[0], handle[1])
    #v2_m = get_measurement_vectors_barnacle(chain_loop, handle[2], handle[3])

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
    adjacent_atoms['P'] = ["O5'", 'OP1', 'OP2']
    adjacent_atoms["O5'"] = ["C5'"]
    adjacent_atoms["C5'"] = ["C4'"]
    adjacent_atoms["C4'"] = ["O4'", "C3'"]
    adjacent_atoms["O4'"] = ["C1'"]
    adjacent_atoms["C1'"] = ["C2'"]
    adjacent_atoms["C2'"] = ["C3'", "O2'"]
    adjacent_atoms["C3'"] = ["O3'"]

    distances = []
    ress = list(chain.get_list())
    for i in range(start_res, end_res+1):
        res = chain[i]
        for key in adjacent_atoms.keys():
            for value in adjacent_atoms[key]:
                distances += [res[key] - res[value]]

    for i in range(start_res, end_res+1):
        distances += [ress[i]['P'] - ress[i-1]["O3'"]]

    return distances


def get_adjacent_interatom_names(chain, start_res, end_res):
    adjacent_atoms = dict()
    adjacent_atoms['P'] = ["O5'", 'OP1', 'OP2']
    adjacent_atoms["O5'"] = ["C5'"]
    adjacent_atoms["C5'"] = ["C4'"]
    adjacent_atoms["C4'"] = ["O4'", "C3'"]
    adjacent_atoms["O4'"] = ["C1'"]
    adjacent_atoms["C1'"] = ["C2'"]
    adjacent_atoms["C2'"] = ["C3'", "O2'"]
    adjacent_atoms["C3'"] = ["O3'"]

    distances = []
    ress = list(chain.get_list())
    for i in range(start_res, end_res+1):
        for key in adjacent_atoms.keys():
            for value in adjacent_atoms[key]:
                distances += [str(key) + "-" + str(value)]

    for i in range(start_res, end_res+1):
        distances += ['%dP-O3' % (i)]

    return distances

def add_residue_to_rosetta_chain(chain, residue):
    '''
    Add a residue and rename all of it's atoms to the Rosetta convention.

    C1' -> C1'

    @param chain: The chain to add to
    @param residue: The residue to be added
    '''
    removed_atoms = []

    for atom in residue.get_list():
        # H-atoms are unnecessary at this moment
        if atom.name.find('H') < 0:
            removed_atoms += [atom]

        residue.detach_child(atom.id)

        atom.name = atom.name.replace("\'", "'")
        atom.id = atom.name

    for atom in removed_atoms:
        residue.add(atom)

    detached_residues = []
    if residue.id[1] in chain:
        detached_residues += [chain[residue.id[1]]]
        #chain.detach_child(chain[residue.id[1]].id)
        chain.detach_child((' ', residue.id[1], ' '))

    chain.add(residue)

    # there should only be one element in the
    # detached residues array
    return detached_residues

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

    fud.pv('handles')
    if handles[1] != length:
        r2_id = chain[handles[1]].id
        chain.detach_child(r2_id)
        loop_chain[handles[3]].id = (' ', handles[1], ' ')
        add_residue_to_rosetta_chain(chain, loop_chain[handles[3]])
    else:
        loop_chain[handles[3]].id = (' ', handles[1], ' ')
        add_residue_to_rosetta_chain(chain, loop_chain[handles[3]])

    # We won't replace the last residue
    # ... or will we?
    counter = 1
    for i in range(handles[2]+1, handles[3]+1):
        loop_chain[i].id = (' ', handles[0] + counter, ' ')
        add_residue_to_rosetta_chain(chain, loop_chain[i])
        counter += 1

def get_initial_measurement_distance(chain_stems, chain_loop, handles):
    '''
    Calculate the rmsd between the measurement vectors after aligning
    the starts.

    @param chain_stems: The PDB coordinates of the chain with the stem.
    @param chain_loop: The PDB coordinates of the sampled loop.
    @param iterations: The number of iterations to use for the CCD loop closure.
    '''
    align_starts(chain_stems, chain_loop, handles)
    c1_target = []
    c1_sampled = []

    for i in range(2):
        target = np.array(get_measurement_vectors(chain_stems, handles[0], handles[1]))
        sampled = np.array(get_measurement_vectors(chain_loop, handles[2], handles[3]))
        '''
        for a in backbone_atoms:
            c1_target += [cuv.magnitude(
                chain_stems[handles[0] - i][a] - chain_stems[handles[1] + i][a])]
            c1_sampled += [cuv.magnitude(
                chain_loop[handles[2] - i][a] - chain_loop[handles[3] + i][a])]
    
    c1_target = np.array(c1_target)
    c1_sampled = np.array(c1_sampled)
        '''
    #return cbc.calc_rmsd(np.array(c1_target), np.array(c1_sampled))
    #return math.sqrt(sum([c ** 2 for c in c1_sampled - c1_target]))
    distances = [cuv.magnitude((sampled - target)[i]) for i in range(len(sampled))]
    rmsd = cuv.vector_set_rmsd(sampled, target)
    #dist = math.sqrt(sum([cuv.magnitude((sampled - target)[i]) ** 2 for i in range(len(sampled))]))
    return rmsd
    #return cuv.magnitude(sampled - target)

def close_fragment_loop(chain_stems, chain_loop, handles, iterations=5000, move_all_angles=True, move_front_angle=True, no_close=False):
    '''
    Align the chain_loop so that it stretches from the end of one stem to the 
    start of the other.

    @param chain_stems: The PDB coordinates of the chain with the stem.
    @param chain_loop: The PDB coordinates of the sampled loop.
    @param iterations: The number of iterations to use for the CCD loop closure.
    '''

    #align_starts(chain_stems, chain_loop, handles)
    e = np.eye(3,3)

    for handle in handles:
        (moving, indeces) = get_atom_coord_array(chain_loop, handle[2], handle[3])
        fixed = np.array(get_measurement_vectors(chain_stems, handle[0], handle[1]))

        start_res = handle[2]
        end_res = handle[3]

        #start_index = indeces[handle[2]+1]
        end_index = indeces[handle[3]+1]

        if no_close:
            rmsd = cbc.calc_rmsd(moving[end_index-3:end_index], fixed)
            return rmsd, chain_loop

        points = []
        #points += [indeces[handle[2]+1]]

        #points += indeces[handle[2]+1] #O3' -> P bond
        if move_all_angles:
            angle_to_move = range(handle[2]+1, handle[3]+1)
        else:
            angle_to_move = [handle[2]+1, handle[3]]

        for i in angle_to_move:
            si = indeces[i]

            # 
            if move_front_angle:
                points += [si]
            points += [si+4, si+5, si+6]

        rot_mat = np.eye(3,3)
            
        '''
        distances = get_adjacent_interatom_distances(chain_loop, handle[2], handle[3])
        names = get_adjacent_interatom_names(chain_loop, handle[2], handle[3])
        '''
        moving = np.array(moving)
        points = np.array(points)

        import fess.aux.ccd.cytvec as cv

        cv.ccd_cython(moving, fixed, points, end_index-3, iterations) 
        rmsd = cbc.calc_rmsd(moving[end_index-3:end_index], fixed)


        chain_loop = set_atom_coord_array(chain_loop, moving, handle[2], handle[3])
        '''
        assert(not np.allclose(moving_orig, moving))

        (moving_new, indeces) = get_atom_coord_array(chain_loop, handles[2], handles[3])

        assert(not np.allclose(moving_orig, moving_new))

        distances2 = get_adjacent_interatom_distances(chain_loop, handles[2], handles[3])

        assert(np.allclose(distances, distances2))
        '''

    return (rmsd, chain_loop)

def align_and_close_loop(seq_len, chain, chain_loop, handles, move_all_angles=True, move_front_angle=True, no_close=False):
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
    
    if handles[0][0] == 0:
        align_starts(chain, chain_loop, handles[0], end=1)
        loop_chain = chain_loop
        r = 0.000
    elif handles[0][1] == seq_len:
        align_starts(chain, chain_loop, handles[0], end=0)
        loop_chain = chain_loop
        r = 0.000
    else:
        r, loop_chain = close_fragment_loop(chain, chain_loop, handles, iterations=10000, move_all_angles=move_all_angles, move_front_angle=move_front_angle, no_close=no_close)

    return (r, loop_chain)

def build_loop(stem_chain, loop_seq, (a,b,i1,i2), seq_len, iterations, consider_contacts=False, consider_starting_pos = True):
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
    import fess.aux.CPDB.BarnacleCPDB as barn

    if consider_contacts:
        model = barn.BarnacleCPDB(loop_seq, 1.9)
    else:
        model = barn.BarnacleCPDB(loop_seq, 0.)

    best_loop_chain = None
    min_energy = (1000000., 100000.)
    prev_energy = min_energy
    handles = (a,b,i1,i2)

    for i in range(iterations):
        sample_len = ss.poisson.rvs(2)
        while sample_len > (len(loop_seq) - 1):
            sample_len = ss.poisson.rvs(2)

        start = rand.randint(0, len(loop_seq) - sample_len)
        end = rand.randint(start + sample_len, len(loop_seq))

        model.sample()
        chain_loop = list(model.structure.get_chains())[0]
        chain_unclosed_loop = copy.deepcopy(chain_loop)

        if handles[0] != 0 and handles[1] != seq_len:
            align_starts(stem_chain, chain_unclosed_loop, [(a,b,i1,i2)], end=2)
        
        loop_chain = copy.deepcopy(chain_unclosed_loop)
        (r, loop_chain) = align_and_close_loop(seq_len, stem_chain, loop_chain, [(a, b, i1, i2)], no_close=False)

        if handles[0] == 0 or handles[1] == seq_len:
            r_start = 0.
        else:
            '''
            r_start = cuv.magnitude(loop_chain[handles[3]]['P'] - 
                                     chain_unclosed_loop[handles[3]]['P'])
            '''
            r_start = cuv.magnitude(loop_chain[handles[2]]['P'] - 
                                     chain_unclosed_loop[handles[3]]['P'])
            

        orig_loop_chain = copy.deepcopy(loop_chain)

        all_chain = copy.deepcopy(stem_chain)
        ftup.trim_chain(loop_chain, i1, i2+1)
        add_loop_chain(all_chain, loop_chain, (a,b,i1,i2), seq_len)

        if consider_contacts:
            contacts2 = ftup.num_noncovalent_clashes(all_chain)
        else:
            contacts2 = 0.

        sys.stderr.write('.')
        sys.stderr.flush()

        if consider_starting_pos:
            energy = (contacts2, r_start)
        else:
            energy = (contacts2, r)

        if energy > prev_energy:
            model.undo()

        prev_energy = energy
        if energy < min_energy:
            min_energy = energy
            min_r = r
            best_loop_chain = copy.deepcopy(orig_loop_chain)
            #output_chain(chain_unclosed_loop, os.path.join(conf.Configuration.test_output_dir, 's3.pdb'))
        '''
        if min_contacts < (0, .1):
            break
        '''

        #trim_chain(loop_chain, i1, i2)

    sys.stderr.write(str(min_energy))
    return (best_loop_chain, min_r)

def reconstruct_loop(chain, sm, ld, side=0, samples=40, consider_contacts=True, consider_starting_pos=True):
    '''
    Reconstruct a particular loop.

    The chain should already have the stems reconstructed.

    @param chain: A Bio.PDB.Chain structure.
    @param sm: A SpatialModel structure
    @param ld: The name of the loop
    '''
    #samples = 2
    bg = sm.bg
    seq = bg.get_flanking_sequence(ld, side)
    (a,b,i1,i2) = bg.get_flanking_handles(ld, side)
    if a == 0 and b == 1:
        # the loop is just a placeholder and doesn't
        # have a length.

        # This should be taken care of in a more elegant
        # manner, but it'll have to wait until it causes
        # a problem
        return None

    # get some diagnostic information
    bl = abs(bg.defines[ld][side * 2 + 1] - bg.defines[ld][side * 2 + 0])
    dist = cuv.vec_distance(bg.coords[ld][1], bg.coords[ld][0])
    dist2 = get_flanking_stem_vres_distance(bg, ld)
    #dist3 = ftug.junction_virtual_atom_distance(bg, ld)
    dist3 = 0.

    sys.stderr.write("reconstructing %s ([%d], %d, %.2f, %.2f, %.2f):" % (ld, len(bg.edges[ld]), bl, dist, dist2, dist3))

    (best_loop_chain, min_r) = build_loop(chain, seq, (a,b,i1,i2), bg.seq_length, samples, consider_contacts, consider_starting_pos)

    #output_chain(chain, os.path.join(conf.Configuration.test_output_dir, 's1.pdb'))
    #output_chain(best_loop_chain, os.path.join(conf.Configuration.test_output_dir, 's2.pdb'))
    print_alignment_pymol_file((a,b,i1,i2))

    ftup.trim_chain(best_loop_chain, i1, i2+1)
    sys.stderr.write('\n')

    add_loop_chain(chain, best_loop_chain, (a,b,i1,i2), bg.seq_length)

    return ((a,b,i1,i2), best_loop_chain, min_r)
    

    #sys.stderr.flush()

from multiprocessing import Process, Pipe
from itertools import izip

def spawn(f):
    def fun(pipe,x):
        pipe.send(f(*x))
        pipe.close()
    return fun

def parmap(f,X):
    pipe=[Pipe() for x in X]
    proc=[Process(target=spawn(f),args=(c,x)) for x,(p,c) in izip(X,pipe)]
    [p.start() for p in proc]
    [p.join() for p in proc]
    return [p.recv() for (p,c) in pipe]

def reconstruct_loops(chain, sm, samples=40, consider_contacts=False):
    '''
    Reconstruct the loops of a chain.

    All of the stems should already be reconstructed in chain.

    @param chain: A Bio.PDB.Chain chain.
    @param sm: The SpatialModel from which to reconstruct the loops.
    '''
    args = []
    for d in sm.bg.defines.keys():
        if d[0] != 's':
            if sm.bg.weights[d] == 2:
                args += [(chain, sm, d, 0, samples, consider_contacts)]
                args += [(chain, sm, d, 1, samples, consider_contacts)]
            else:
                args += [(chain, sm, d, 0, samples, consider_contacts)]

    #pool = mp.Pool(processes=4)
    #r = parmap(reconstruct_loop, args)
    r = [reconstruct_loop(*arg) for arg in args]
    r = [x for x in r if x is not None]

    for ((a,b,i1,i2), best_loop_chain, min_r) in r:
        add_loop_chain(chain, best_loop_chain, (a,b,i1,i2), sm.bg.length)

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
    av = { 'U': ['N1', "C1'", "C2'"], 'C': ['N1', "C1'", "C2'"], 'A': ['N9', "C1'", "C2'"], 'G': ['N9', "C1'", "C2'"], 'rU': ['N1', "C1'", "C2'"], 'rC': ['N1', "C1'", "C2'"], 'rA': ['N9', "C1'", "C2'"], 'rG': ['N9', "C1'", "C2'"] }

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

    ress = chain.get_list()

    for i in range(len(ress)):
        num = ress[i].id[1]
        name = ress[i].resname.strip()

        if num-1 == len(seq):
            continue

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
def mend_breakpoint(h, chain, source_chain):
    # try to make the last residue of the connecting stem transition
    # into the loop region

    # this is done by making a copy of the loop region and aligning
    # it to the backbone of the last residue of the stem

    # the new section (last residue of the stem and two residues of
    # the loop, newly aligned) are then loop-closed to align to the
    # original orientation of the fitted loop region
    temp_loop_chain = copy.deepcopy(source_chain)
    align_starts(chain, temp_loop_chain, [h], end=2)
    rev_handles = [(h[2]-1, h[2]+1, h[0]-1, h[0]+1)]
    temp_loop_chain[h[2] + 1].id = (' ', h[0]+1, ' ')
    temp_loop_chain[h[2] + 2].id = (' ', h[0]+2, ' ')
    temp_loop_chain[h[2] + 3].id = (' ', h[0]+3, ' ')

    detached_residues = []
    detached_residues += add_residue_to_rosetta_chain(chain, temp_loop_chain[h[2]+1])
    detached_residues += add_residue_to_rosetta_chain(chain, temp_loop_chain[h[2]+2])
    #detached_residues += add_residue_to_rosetta_chain(chain, temp_loop_chain[h[2]+3])
    (r1, stem_chain) = align_and_close_loop(10000, source_chain, chain, rev_handles, move_all_angles=False, move_front_angle=False)

    for dr in detached_residues:
        add_residue_to_rosetta_chain(chain, dr)

def get_ideal_stem(template_stem_length):
    '''
    Load the model for an ideal stem of a particular length.

    TODO: Remove the duplicate code in graph_pdb.py (get_mids_core)
    '''
    tstart1 = 1
    tstart2 = template_stem_length * 2
    tend1 = template_stem_length
    tend2 = template_stem_length + 1

    template_filename = 'ideal_1_%d_%d_%d.pdb' % (tend1, tend2, tstart2)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ideal_chain = list(bpdb.PDBParser().get_structure('test', 
                op.join(conf.Configuration.stem_fragment_dir, template_filename)).get_chains())[0]

        chain = ftug.extract_define_residues([tstart1,tend1,tend2,tstart2], ideal_chain)
    return chain

def mend_breakpoint_new(chain, res1, res2):
    '''
    Try and mend the breakpoint between res1 and res2.

    This will be done by excising all residues in the range [res1, res2].
    The excised region will be replaced by a new, connected region (nr),
    consisting of residues [res1-1, res2+1]. The new region will be aligned
    along residue res1-1 and will be loop-closed to res2+1. It will
    then be inserted into the chain as the new section.
    '''
    region_len = res2 - res1 + 2
    replacement_length = get_ideal_stem(region_len)
    h = [(res1,res2+1,1,region_len)]

    nr = get_ideal_stem(region_len)
    align_starts(chain, nr, h, end=2)

    h1 = h

    (r1, loop_chain) = align_and_close_loop(10000, chain, nr, h, move_all_angles=False, move_front_angle=False)
    add_loop_chain(chain, nr, h[0], h[0][3] - h[0][2])

def align_source_to_target_fragment(target_chain, source_chain, sm, angle_def, ld):
    '''
    Align a PDB chain to the position where it is supposed to
    bridge the gap between two stems.

    @param target_chain: The chain the fragment is to be inserted into
    @param source_chain: The chain from which the fragment comes
    @param sm: The SpatialModel being reconstructed
    @param angle_def: The define containing the residue numbers in source_chain
    @param ld: The name of the fragment.
    '''
    connections = sm.bg.connections(ld)

    (s1b, s1e) = sm.bg.get_sides(connections[0], ld)
    #(s2b, s2e) = sm.bg.get_sides(connections[1], ld)

    (sd, bd) = sm.bg.get_sides_plus(connections[0], ld)

    t_v = (target_chain[sm.bg.defines[connections[0]][sd]]["C3'"].get_vector().get_array(),
           target_chain[sm.bg.defines[connections[0]][sd]]["C4'"].get_vector().get_array(), 
           target_chain[sm.bg.defines[connections[0]][sd]]["O4'"].get_vector().get_array())

    s_v = (source_chain[angle_def.define[bd]]["C3'"].get_vector().get_array(),
           source_chain[angle_def.define[bd]]["C4'"].get_vector().get_array(),
           source_chain[angle_def.define[bd]]["O4'"].get_vector().get_array())

    t_centroid = cuv.get_vector_centroid(t_v)
    s_centroid = cuv.get_vector_centroid(s_v)

    t_v1 = t_v - t_centroid
    s_v1 = s_v - s_centroid

    sup = brmsd.optimal_superposition(s_v1, t_v1)

    for atom in bpdb.Selection.unfold_entities(source_chain, 'A'):
        atom.transform(np.eye(3,3), -s_centroid)
        atom.transform(sup, t_centroid)

    pass

def reconstruct_bulge_with_fragment_core(chain, source_chain, sm, ld, sd0, sd1, angle_def, move_all_angles=False):
    connections = sm.bg.connections(ld)

    a0 = sm.bg.defines[connections[0]][sd0]
    b0 = sm.bg.defines[connections[1]][sd1]
    a = [a0,b0]
    a.sort()
    (a0,b0) = a

    if len(angle_def.define) == 4:
        '''
        a0_1 = sm.bg.defines[connections[0]][sm.bg.same_stem_end(sd0)]
        b0_1 = sm.bg.defines[connections[1]][sm.bg.same_stem_end(sd1)]
        b = [a0_1, b0_1]
        b.sort()
        (a0_1, b0_1) = b
        '''

        # sort the defines by the first entry in each define
        # i.e. [3,4,1,2] -> [1,2,3,]
        s1 = map(list, zip(*[iter(sm.bg.defines[ld])]*2))
        s1.sort()

        for s in s1:
            s[0] -= 1
            s[1] += 1

        s2 = map(list, zip(*[iter(angle_def.define)]*2))
        s2.sort()

        for s in s2:
            s[0] -= 1
            s[1] += 1

        # Associate the defines of the source with those of the loop
        # according to which are lower and which are greater:
        # s1: [(21, 22), (46, 48)]
        # s2: [(134, 135), (142, 144)]
        # handles: [[21, 22, 134, 135], [46, 48, 142, 144]]

        handles = [[i for l in s for i in l] for s in zip(s1, s2)]
        '''

        if angle_def.ang_type == 2 and sd0 == 1:
        #if False:
            i1_0 = angle_def.define[2]
            i2_0 = angle_def.define[3]
            i1_1 = angle_def.define[0]
            i2_1 = angle_def.define[1]
        else:
            i1_0 = angle_def.define[0]
            i2_0 = angle_def.define[1]
            i1_1 = angle_def.define[2]
            i2_1 = angle_def.define[3]

        handles = [(a0,b0,i1_0,i2_0), (a0_1, b0_1, i1_1, i2_1)]
        '''
    else:
        i1_0 = angle_def.define[0]
        i2_0 = angle_def.define[1]
        handles = [(a0,b0,i1_0-1,i2_0+1)]

    seq_len = handles[0][3] - handles[0][2] # i2_0 - i1_0
    align_starts(chain, source_chain, handles, end=0)

    (r, loop_chain) = align_and_close_loop(seq_len, chain, source_chain, handles, move_all_angles=move_all_angles)

    for h in handles:
        mend_breakpoint(h, chain, source_chain)
    for h in handles:
        add_loop_chain(chain, source_chain, h, h[3] - h[2])
        mend_breakpoint_new(chain, h[0], h[0]+1)

    #(r, chain) = align_and_close_loop(seq_len, source_chain, chain, 

def reconstruct_bulge_with_fragment(chain, sm, ld, fragment_library=dict(), move_all_angles=False):
    '''
    Reconstruct a loop with the fragment its statistics were derived from.

    @param chain: The chain containing the reconstructed stems.
    @param sm: The SpatialModel containing the information about the sampled
        stems and angles
    @param ld: The name of the loop to reconstruct.
    '''

    #find some potential sides
    #both ways should work
    #i.e. if [0][1] is present, [0][1] should also be present
    for key1 in sm.angle_defs[ld].keys():
        break

    angle_def = sm.angle_defs[ld][key1]

    # the file containing the pdb coordinates of this fragment
    filename = '%s_%s.pdb' % (angle_def.pdb_name, "_".join(map(str, angle_def.define)))
    filename = os.path.join(conf.Configuration.stem_fragment_dir, filename)

    # do some caching while loading the filename
    if filename in fragment_library.keys():
        source_chain = fragment_library[filename].copy()
    else:
        source_chain = list(bpdb.PDBParser().get_structure('temp', filename).get_chains())[0]
        fragment_library[filename] = source_chain

    #align_source_to_target_fragment(chain, source_chain, sm, angle_def, ld)

    connections = sm.bg.connections(ld)
    (sd0, bd0) = sm.bg.get_sides_plus(connections[0], ld)
    (sd1, bd1) = sm.bg.get_sides_plus(connections[1], ld)


    if len(angle_def.define) == 2:
        reconstruct_bulge_with_fragment_core(chain, source_chain, sm, ld, sd0, sd1, angle_def, move_all_angles=move_all_angles)
        return

    reconstruct_bulge_with_fragment_core(chain, source_chain, sm, ld, sd0, sd1, angle_def, move_all_angles=move_all_angles)

    return

def reconstruct_loop_with_fragment(chain, sm, ld, fragment_library=dict()):
    '''
    Reconstruct a loop with the fragment its statistics were derived from.

    @param chain: The chain containing the reconstructed stems.
    @param sm: The SpatialModel containing the information about the sampled
        stems and angles
    @param ld: The name of the loop to reconstruct.
    '''

    loop_def = sm.loop_defs[ld]
    angle_def = loop_def

    fud.pv('ld')
    fud.pv('loop_def')

    if loop_def.define[1] - loop_def.define[0] == 1:
        return

    # the file containing the pdb coordinates of this fragment
    filename = '%s_%s.pdb' % (loop_def.pdb_name, "_".join(map(str, loop_def.define)))
    filename = os.path.join(conf.Configuration.stem_fragment_dir, filename)

    # do some caching while loading the filename
    if filename in fragment_library.keys():
        source_chain = fragment_library[filename].copy()
    else:
        source_chain = list(bpdb.PDBParser().get_structure('temp', filename).get_chains())[0]
        fragment_library[filename] = source_chain

    align_source_to_target_fragment(chain, source_chain, sm, loop_def, ld)
    connection = list(sm.bg.edges[ld])[0]

    (sd0, bd0) = sm.bg.get_sides_plus(connection, ld)

    fud.pv('sm.bg.defines[connection]')
    if sd0 == 0:
        a0,b0 = sm.bg.defines[connection][0], sm.bg.defines[connection][3]
    else:
        a0,b0 = sm.bg.defines[connection][1], sm.bg.defines[connection][2]

    i1_0 = angle_def.define[0] - 1
    i2_0 = angle_def.define[1] + 1

    seq_len = i2_0 - i1_0
    align_starts(chain, source_chain, [(a0,b0,i1_0,i2_0)], end=0)
    (r, loop_chain) = align_and_close_loop(seq_len, chain, source_chain, [(a0,b0,i1_0,i2_0)], move_all_angles=False)
    add_loop_chain(chain, source_chain, (a0,b0,i1_0,i2_0), i2_0 - i1_0)
    mend_breakpoint_new(chain, a0, a0+1)

    return

def reconstruct_fiveprime_with_fragment(chain, sm, ld, fragment_library=dict()):
    '''
    Reconstruct a 5' unpaired region with the fragment its statistics were derived from.

    @param chain: The chain containing the reconstructed stems.
    @param sm: The SpatialModel containing the information about the sampled
        stems and angles
    @param ld: The name of the loop to reconstruct.
    '''

    try:
        fiveprime_def = sm.fiveprime_defs[ld]
    except:
        reconstruct_loop(chain, sm, ld)
        return

    if fiveprime_def.define[1] - fiveprime_def.define[0] == 1:
        return


    # the file containing the pdb coordinates of this fragment
    filename = '%s_%s.pdb' % (fiveprime_def.pdb_name, "_".join(map(str, fiveprime_def.define)))
    filename = os.path.join(conf.Configuration.stem_fragment_dir, filename)

    # do some caching while loading the filename
    if filename in fragment_library.keys():
        source_chain = fragment_library[filename].copy()
    else:
        source_chain = list(bpdb.PDBParser().get_structure('temp', filename).get_chains())[0]
        fragment_library[filename] = source_chain

    align_source_to_target_fragment(chain, source_chain, sm, fiveprime_def, ld)

    # add the new chain to the old one
    for j in range(0, len(fiveprime_def.define), 2):
        for k in range(max(fiveprime_def.define[j],1), fiveprime_def.define[j+1]+1):
            target_index = sm.bg.defines[ld][j] + k - fiveprime_def.define[j]

            if target_index in chain:
                chain.detach_child(chain[target_index].id)

            e = source_chain[k]
            e.id = (e.id[0], target_index, e.id[2])

            chain.add(e)
    pass

def reconstruct_threeprime_with_fragment(chain, sm, ld, fragment_library=dict()):
    '''
    Reconstruct a 5' unpaired region with the fragment its statistics were derived from.

    @param chain: The chain containing the reconstructed stems.
    @param sm: The SpatialModel containing the information about the sampled
        stems and angles
    @param ld: The name of the loop to reconstruct.
    '''

    threeprime_def = sm.threeprime_defs[ld]

    # the file containing the pdb coordinates of this fragment
    filename = '%s_%s.pdb' % (threeprime_def.pdb_name, "_".join(map(str, threeprime_def.define)))
    filename = os.path.join(conf.Configuration.stem_fragment_dir, filename)

    # do some caching while loading the filename
    if filename in fragment_library.keys():
        source_chain = fragment_library[filename].copy()
    else:
        source_chain = list(bpdb.PDBParser().get_structure('temp', filename).get_chains())[0]
        fragment_library[filename] = source_chain

    align_source_to_target_fragment(chain, source_chain, sm, threeprime_def, ld)

    # add the new chain to the old one
    for j in range(0, len(threeprime_def.define), 2):
        for k in range(threeprime_def.define[j], threeprime_def.define[j+1]+1):
            target_index = sm.bg.defines[ld][j] + k - threeprime_def.define[j]

            if target_index in chain:
                chain.detach_child(chain[target_index].id)

            e = source_chain[k]
            e.id = (e.id[0], target_index, e.id[2])

            chain.add(e)
    pass

