import multiprocessing as mp

import corgy.builder.models as models
import corgy.builder.rmsd as brmsd

import corgy.graph.graph_pdb as cgg
import corgy.utilities.vector as cuv
import corgy.utilities.debug as cud
import corgy.utilities.pdb as cup

import corgy.builder.ccd as cbc
import corgy.aux.ccd.cytvec as cv

#import corgy.aux.Barnacle as barn
import corgy.aux.CPDB.src.examples.BarnacleCPDB as barn

import Bio.PDB as bpdb
import Bio.PDB.Chain as bpdbc
import Bio.PDB.Model as bpdbm
import Bio.PDB.Structure as bpdbs

import scipy.stats as ss

from scipy.stats import norm, poisson

import os, math, sys
import corgy.builder.config as conf
import copy, time
import random as rand

import numpy as np

def get_measurement_vectors1(ress, r1, r2):
    return( ress[r2]['C4*'].get_vector().get_array(), 
            ress[r2]['C3*'].get_vector().get_array(),
            ress[r2]['O3*'].get_vector().get_array())

def get_measurement_vectors2(ress, r1, r2):
    return( ress[r2]['O4*'].get_vector().get_array(), 
            ress[r2]['C1*'].get_vector().get_array(),
            ress[r2]['C2*'].get_vector().get_array())

def pdb_rmsd(c1, c2, backbone=True, superimpose=True):
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
        cud.pv('len(c1.get_list())')
        cud.pv('len(c2.get_list())')
        print >>sys.stderr, "Chains of different length"
        raise Exception("Chains of different length.")

    #for i in range(1, len(list(c1.get_list()))+1):
    for i in [r.id[1] for r in c1.get_residues()]:
        if backbone:
            anames = a_5_names + a_names[c1[i].resname.strip()] + a_3_names
        else:
            anames = a_5_names + a_3_names
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

    #print "rmsd len:", len(all_atoms1), len(all_atoms2)
    if superimpose:
        sup = bpdb.Superimposer()
        sup.set_atoms(all_atoms1, all_atoms2)

        sup.apply(c2.get_atoms())

        return (len(all_atoms1), sup.rms)
    else:
        crvs1 = np.array([a.get_vector().get_array() for a in all_atoms1])
        crvs2 = np.array([a.get_vector().get_array() for a in all_atoms2])

        return (len(all_atoms1), brmsd.rmsd(crvs1, crvs2))

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
            (vr1_p, vr1_v, vr1_v_l, vr1_v_r) = cgg.virtual_res_3d_pos(bg, connecting_stems[0], bg.stem_length(connecting_stems[0]) - 1)
        else:
            (vr1_p, vr1_v, vr1_v_l, vr1_v_r) = cgg.virtual_res_3d_pos(bg, connecting_stems[0], 0)

        if s2b == 1:
            (vr2_p, vr2_v, vr2_v_l, vr2_v_r) = cgg.virtual_res_3d_pos(bg, connecting_stems[1], bg.stem_length(connecting_stems[1]) - 1)
        else:
            (vr2_p, vr2_v, vr2_v_l, vr2_v_r) = cgg.virtual_res_3d_pos(bg, connecting_stems[1], 0)

        dist2 = cuv.vec_distance((vr1_p + 7 * vr1_v), (vr2_p + 7. * vr2_v))
    else:
        dist2 = 0.

    return dist2

a_5_names = ['P', 'OP1', 'OP2', 'P', 'O5*', 'C5*', 'C4*', 'O4*', 'O2*']
a_3_names = ['C1*', 'C2*', 'C3*', 'O3*']

backbone_atoms = ['P', 'O5*', 'C5*', 'C4*', 'C3*', 'O3*']

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

    #start_index = indeces[handles[2]+1]
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
    if consider_contacts:
        model = barn.BarnacleCPDB(loop_seq, 1.9)
    else:
        model = barn.BarnacleCPDB(loop_seq, 0.)

    best_loop_chain = None
    min_energy = (1000000., 100000.)
    prev_energy = min_energy
    handles = (a,b,i1,i2)

    cud.pv('cup.num_noncovalent_clashes(stem_chain)')

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
            align_starts(stem_chain, chain_unclosed_loop, (a,b,i1,i2), end=0)
        
        (r, loop_chain) = align_and_close_loop(seq_len, stem_chain, chain_loop, (a, b, i1, i2))
        cud.pv('r')
        if handles[0] == 0 or handles[1] == seq_len:
            r_start = 0.
        else:
            r_start = cuv.magnitude(loop_chain[handles[2]+2]['P'] - 
                                     chain_unclosed_loop[handles[2]+2]['P'])

        orig_loop_chain = copy.deepcopy(loop_chain)

        all_chain = copy.deepcopy(stem_chain)
        cup.trim_chain(loop_chain, i1, i2+1)
        add_loop_chain(all_chain, loop_chain, (a,b,i1,i2), seq_len)

        if consider_contacts:
            contacts2 = cup.num_noncovalent_clashes(all_chain)
        else:
            contacts2 = 0.

        sys.stderr.write('.')
        sys.stderr.flush()

        energy = (contacts2, r_start * r)
        cud.pv('(start,end,energy)')
        if energy > prev_energy:
            model.undo()

        prev_energy = energy
        if energy < min_energy:
            min_energy = energy
            min_r = r
            best_loop_chain = copy.deepcopy(orig_loop_chain)
            output_chain(chain_unclosed_loop, os.path.join(conf.Configuration.test_output_dir, 's3.pdb'))
            cud.pv('min_energy')
        '''
        if min_contacts < (0, .1):
            break
        '''

        #trim_chain(loop_chain, i1, i2)

    sys.stderr.write(str(min_energy))
    return (best_loop_chain, min_r)

def reconstruct_loop(chain, sm, ld, side=0, samples=40, consider_contacts=True):
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

    # get some diagnostic information
    bl = abs(bg.defines[ld][side * 2 + 1] - bg.defines[ld][side * 2 + 0])
    dist = cuv.vec_distance(bg.coords[ld][1], bg.coords[ld][0])
    dist2 = get_flanking_stem_vres_distance(bg, ld)

    sys.stderr.write("reconstructing %s ([%d], %d, %f, %f):" % (ld, len(bg.edges[ld]), bl, dist, dist2))
    '''
    if a == 0 and b == 1:
        # the loop is just a placeholder and doesn't
        # have a length.

        # This should be taken care of in a more elegant
        # manner, but it'll have to wait until it causes
        # a problem
        return
    '''

    (best_loop_chain, min_r) = build_loop(chain, seq, (a,b,i1,i2), bg.length, samples, consider_contacts)

    output_chain(chain, os.path.join(conf.Configuration.test_output_dir, 's1.pdb'))
    output_chain(best_loop_chain, os.path.join(conf.Configuration.test_output_dir, 's2.pdb'))
    print_alignment_pymol_file((a,b,i1,i2))

    cup.trim_chain(best_loop_chain, i1, i2+1)
    return ((a,b,i1,i2), best_loop_chain, min_r)
    
    add_loop_chain(chain, best_loop_chain, (a,b,i1,i2), bg.length)

    sys.stderr.write('\n')
    sys.stderr.flush()

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
    r = parmap(reconstruct_loop, args)
    #r = [reconstruct_loop(*arg) for arg in args]
    for ((a,b,i1,i2), best_loop_chain, min_r) in r:
        add_loop_chain(chain, best_loop_chain, (a,b,i1,i2), sm.bg.length)

    cud.pv('r')

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

    t_v = (target_chain[sm.bg.defines[connections[0]][sd]]['C3*'].get_vector().get_array(),
           target_chain[sm.bg.defines[connections[0]][sd]]['C4*'].get_vector().get_array(), 
           target_chain[sm.bg.defines[connections[0]][sd]]['O4*'].get_vector().get_array())

    s_v = (source_chain[angle_def.define[bd]]['C3*'].get_vector().get_array(),
           source_chain[angle_def.define[bd]]['C4*'].get_vector().get_array(),
           source_chain[angle_def.define[bd]]['O4*'].get_vector().get_array())

    t_centroid = cuv.get_vector_centroid(t_v)
    s_centroid = cuv.get_vector_centroid(s_v)

    t_v1 = t_v - t_centroid
    s_v1 = s_v - s_centroid

    sup = brmsd.optimal_superposition(s_v1, t_v1)

    for atom in bpdb.Selection.unfold_entities(source_chain, 'A'):
        atom.transform(np.eye(3,3), -s_centroid)
        atom.transform(sup, t_centroid)

    pass

def reconstruct_bulge_with_fragment(chain, sm, ld, fragment_library=dict()):
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
        for key2 in sm.angle_defs[ld][key1].keys():
            break

    angle_def = sm.angle_defs[ld][key1][key2]

    # the file containing the pdb coordinates of this fragment
    filename = '%s_%s.pdb' % (angle_def.pdb_name, "_".join(map(str, angle_def.define)))
    filename = os.path.join(conf.Configuration.stem_fragment_dir, filename)

    # do some caching while loading the filename
    if filename in fragment_library.keys():
        source_chain = fragment_library[filename].copy()
    else:
        source_chain = list(bpdb.PDBParser().get_structure('temp', filename).get_chains())[0]
        fragment_library[filename] = source_chain

    align_source_to_target_fragment(chain, source_chain, sm, angle_def, ld)

    connections = sm.bg.connections(ld)
    (s1b, s1e) = sm.bg.get_sides(connections[0], ld)
    (sd, bd) = sm.bg.get_sides_plus(connections[0], ld)

    d = angle_def.define
    d1 = sm.bg.defines[ld]

    # add the new chain to the old one
    for j in range(0, len(angle_def.define), 2):
        #bd = sm.bg.get_bulge_dimensions(ld)
        if len(d) == 4:
            bd = (d1[1] - d1[0], d1[3] - d1[2])
        for k in range(angle_def.define[j], angle_def.define[j+1]+1):
            target_index = sm.bg.defines[ld][j] + k - angle_def.define[j]
            if len(d) == 4:
                if bd != (d[1] - d[0], d[3] - d[2]):
                    target_index = sm.bg.defines[ld][(j+2) % 4] + k - angle_def.define[j]

            if target_index in chain:
                #print >> sys.stderr, "detaching...", target_index
                chain.detach_child(chain[target_index].id)

            e = source_chain[k]
            e.id = (e.id[0], target_index, e.id[2])

            chain.add(e)
    pass

def reconstruct_loop_with_fragment(chain, sm, ld, fragment_library=dict()):
    '''
    Reconstruct a loop with the fragment its statistics were derived from.

    @param chain: The chain containing the reconstructed stems.
    @param sm: The SpatialModel containing the information about the sampled
        stems and angles
    @param ld: The name of the loop to reconstruct.
    '''

    loop_def = sm.loop_defs[ld]

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

    # add the new chain to the old one
    for j in range(0, len(loop_def.define), 2):
        for k in range(loop_def.define[j], loop_def.define[j+1]+1):
            target_index = sm.bg.defines[ld][j] + k - loop_def.define[j]

            if target_index in chain:
                #print >> sys.stderr, "detaching...", target_index
                chain.detach_child(chain[target_index].id)

            e = source_chain[k]
            e.id = (e.id[0], target_index, e.id[2])

            chain.add(e)
    pass

def reconstruct_fiveprime_with_fragment(chain, sm, ld, fragment_library=dict()):
    '''
    Reconstruct a 5' unpaired region with the fragment its statistics were derived from.

    @param chain: The chain containing the reconstructed stems.
    @param sm: The SpatialModel containing the information about the sampled
        stems and angles
    @param ld: The name of the loop to reconstruct.
    '''

    fiveprime_def = sm.fiveprime_defs[ld]

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
                #print >> sys.stderr, "detaching...", target_index
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
                #print >> sys.stderr, "detaching...", target_index
                chain.detach_child(chain[target_index].id)

            e = source_chain[k]
            e.id = (e.id[0], target_index, e.id[2])

            chain.add(e)
    pass
