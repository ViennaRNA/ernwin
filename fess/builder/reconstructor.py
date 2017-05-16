import multiprocessing as mp
import warnings

import itertools as it
import fess.builder.models as models
import fess.builder.config as fbc

import os.path as op

import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.pdb as ftup
import forgi.threedee.utilities.average_atom_positions as ftua
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.vector as cuv
ftuv = cuv
import forgi.utilities.debug as fud

import fess.builder.ccd as cbc

import forgi.threedee.model.similarity as brmsd
import forgi.graph.bulge_graph as fgb
#import borgy.aux.Barnacle as barn
#import fess.aux.CPDB.BarnacleCPDB as barn

import Bio.PDB as bpdb
import Bio.PDB.Atom as bpdba
import Bio.PDB.Residue as bpdbr
import Bio.PDB.Chain as bpdbc
import Bio.PDB.Model as bpdbm
import Bio.PDB.Structure as bpdbs

import scipy.stats as ss

from scipy.stats import norm, poisson

import os, math, sys
import fess.builder.config as conf
import copy, time
import random as rand

import logging
log=logging.getLogger(__name__)

import numpy as np

class Reconstructor(object):
    def __init__(self, pdb_library_path, cg_library_path):
        self.cg_library_path = cg_library_path
        self.pdb_library_path = pdb_library_path
        self.stem_use_average_method = True
    def reconstruct(self, sm):
        '''
        Re-construct a full-atom model from a coarse-grain model.

        :param sm: The Spatial Model to reconstruct
        :returns: A dictionary {chain_id : Bio.PDB.Chain.Chain}
        '''
        # First, reconstruct the stems, which will be used as scaffold for the bulges
        chains = {}
        for stem in sm.elem_defs:
            if stem[0] != "s":
                continue
            self._reconstruct_stem(sm, stem, chains)

        #Some validation
        _compare_cg_chains_partial(sm.bg, chains)
        ftup.output_multiple_chains(chains.values(), "reconstr_stems.pdb")
        '''# Some more validation. Useless unless searching for a bug
        chains_in_file = ftup.get_all_chains("reconstr_stems.pdb")
        chains_in_file = { c.id:c for c in chains_in_file }
        for c in chains:
            r = ftup.pdb_rmsd(chains[c], chains_in_file[c], superimpose = False)[1]
            assert r<0.001, "r={} for chain {}".format(r, c)
        '''

        for loop in sm.bg.defines:
            if loop[0]=="s": continue
            self._reconstruct_with_fragment(chains, sm, loop)
        ftup.output_multiple_chains(chains.values(), "reconstr_all.pdb")

        replace_bases(chains, sm.bg)
        ftup.output_multiple_chains(chains.values(), "base_replaced.pdb")
        reorder_residues(chains, sm.bg)
        return chains

    def _get_source_cg_and_chain(self, stat):
        """
        Load the fragment defined in the stat from the fragment library as pdb and cg.

        :param stat: The forgi.threedee.model.stats.StemStat or ftms.AngleStat or ftms.LoopStat object.
        """
        stat_name = stat.pdb_name
        pdb_basename = stat_name.split(":")[0]
        pdb_filename = op.expanduser(op.join(self.pdb_library_path, "".join(pdb_basename.split("_")[:-1])+".pdb"))
        cg_filename = op.expanduser(op.join(self.cg_library_path, pdb_basename+".cg"))
        #Make sure the files exist.
        with open(pdb_filename): pass
        with open(cg_filename): pass
        log.debug("Opening cg-file %s to extract stat %s", cg_filename, stat.pdb_name)
        cg = ftmc.CoarseGrainRNA(cg_filename) #The cg with the template

        chains = ftup.get_all_chains(pdb_filename)
        chains = {c.id:c for c in chains}
        return cg, chains

    def _reconstruct_stem(self, sm, stem_name, new_chains):
        '''
        Reconstruct a particular stem.

        :param sm: The SpatialModel that should be translated to PDB
        :param stem_name: The name of the stem to be reconstructed.
        :param new_chains: A dict chainid:chain that will be filled with the reconstructed model
        '''
        stem = sm.stems[stem_name]
        stem_stat = sm.elem_defs[stem_name]
        orig_def = sm.bg.defines[stem_name]
        cg_orig = sm.bg
        cg, chains = self._get_source_cg_and_chain(stem_stat)
        #: The define of the stem in the original cg where the fragment is from
        sd = cg.get_node_from_residue_num(stem_stat.define[0])
        chains  = ftup.extract_subchains_from_seq_ids(chains,
                                                      cg.define_residue_num_iterator(sd, seq_ids=True))
        if stem_stat.define != cg.defines[sd]:
            log.error("%s != %s for %s (%s)", stem_stat.define, cg.defines[sd], sd, stem_stat.pdb_name)
            raise ValueError("The CG files where the stats where extracted and the cg file used for reconstruction are not consistent!")

        _align_chain_to_stem(cg, chains, sd, stem, self.stem_use_average_method)

        for i in range(stem_stat.bp_length):
            for strand in range(2):
                target_resid = cg_orig.seq_ids[ orig_def[strand*2] + i - 1 ]
                source_resid = cg.seq_ids[ stem_stat.define[strand*2] + i - 1 ]
                residue = chains[source_resid.chain][source_resid.resid]
                #Change the resid to the target
                residue.id = target_resid.resid
                assert residue.id == target_resid.resid
                if target_resid.chain not in new_chains:
                    log.info("Adding chain with id %r for residue %r", target_resid.chain, target_resid)
                    new_chains[target_resid.chain] =  bpdb.Chain.Chain(target_resid.chain)
                #Now, add the residue to the target chain
                new_chains[target_resid.chain].add(residue)

        assert _validate_pdb_to_stem(_define_to_stem_model(cg, chains, sd), new_chains, cg_orig, cg_orig.get_node_from_residue_num(orig_def[0]))

        return new_chains

    def _reconstruct_with_fragment(self, chains, sm, ld):
        '''
        Reconstruct a loop with the fragment its statistics were derived from.

        :param chains: A dict chain_id:chain that will be filled.
        :param sm: The SpatialModel containing the information about the sampled
            stems and angles
        :param ld: The name of the loop to reconstruct.

        :returns: Modifies chains in-place and returns a reference to it.
        '''

        close_loop = True

        try:
            angle_stat = sm.elem_defs[ld]
        except KeyError:
            # Not part of the minimal spanning tree. We probably need to do use CCD or BARNACLE
            warnings.warn("ML that is not part o the MST is not yet implemented!")
            return

        cg_from, chains_from = self._get_source_cg_and_chain(angle_stat)
        cg_to = sm.bg
        chains_to = chains
        elem_to = ld
        elem_from = cg_from.get_node_from_residue_num(angle_stat.define[0])

        insert_element(cg_to, cg_from, elem_to, elem_from, chains_to, chains_from)

        return chains

def _align_chain_to_stem(cg, chains, elem_name, stem2, use_average_method=True):
    """
    Rotate and tranlate chains to match the orientation of the coarse-grained stem.

    :param cg: The coarse-grained RNA where the fragment is originally from.
    :param chains: The PDB chains containing the original fragment.
    :param elem_name: The element name of the stem in cg.
    :param stem2: The target (cg-)stem. A StemModel object.
    """
    #The stem fragment we will rotate and translate
    stem1 = _define_to_stem_model(cg, chains, elem_name)

    '''
    (r, u, v, t) = ftug.get_stem_orientation_parameters(stem1.vec(),
                                                       (stem1.twists[0] + stem1.twists[1]) / 2.,
                                                       stem2.vec(),
                                                       (stem2.twists[0] + stem2.twists[1]) / 2.)
    '''
    if not use_average_method:
        (r, u, v, t) = ftug.get_stem_orientation_parameters(stem1.vec(),
                                                           stem1.twists[0],
                                                           stem2.vec(),
                                                           stem2.twists[0])
    else:
        tw1 = ftug.virtual_res_3d_pos_core(stem1.mids, stem1.twists, 2, 4)[1]
        tw2 = ftug.virtual_res_3d_pos_core(stem2.mids, stem2.twists, 2, 4)[1]
        (r, u, v, t) = ftug.get_stem_orientation_parameters(stem1.vec(),
                                                           tw1,
                                                           stem2.vec(),
                                                           tw2)
    rot_mat = get_stem_rotation_matrix(stem1, stem2, use_average_method)

    assert np.allclose(ftuv.normalize(stem2.vec()),  ftuv.normalize(np.dot(stem1.vec(), rot_mat)))
    rotate_chain(chains, rot_mat, (stem1.mids[0] + stem1.mids[1]) / 2.)
    translate_chain(chains, (stem2.mids[0] + stem2.mids[1]) / 2. - (stem1.mids[0] + stem1.mids[1]) / 2.)

    assert _validate_pdb_to_stem(stem2, chains, cg, elem_name)

def _define_to_stem_model(cg, chains, elem_name):
    '''
    Extract a StemModel from a Bio.PDB.Chain structure.

    @param chain: The Bio.PDB.Chain representation of the chain
    @param define: The BulgeGraph element name, e.g. "s0"
    @return: A StemModel with the coordinates and orientation of the stem.
    '''
    stem = models.StemModel(name=elem_name)

    coords, twists = ftug.stem_from_chains(cg, chains, elem_name)

    stem.mids = coords
    stem.twists = twists

    # Some validation:
    atom_coords = []
    for chain in chains.values():
        for atom in bpdb.Selection.unfold_entities(chain, 'A'):
            atom_coords.append(atom.coord)
    s = np.sum(atom_coords, axis=0)
    s = s/len(atom_coords)
    log.debug("Center of atoms: %s, Center of mids: %s. Stem mids: %s", s, (stem.mids[1]+stem.mids[0])/2, stem.mids)
    return stem

def _validate_pdb_to_stem(target_stem, chains, cg, elem_name):
    """
    :param target_stem: A StemModel to which the pdb chain should be aligned
    :param chains: A dict {chain_id: Chain}
    :param cg: The original coarse-grained representation of the pdb chains
    :param elem_name: The elem_name in cg.
    """
    try:
        pdb_stem = _define_to_stem_model(cg, chains, elem_name)
    except:
        for chain in chains.values():
            log.error([r.id for r in chain.get_residues()])
        raise
    d_start = ftuv.magnitude(pdb_stem.mids[0] - target_stem.mids[0])
    d_end   = ftuv.magnitude(pdb_stem.mids[1] - target_stem.mids[1])
    assert d_start<0.1, "{Distance between stem starts {} and {} is too big: {}".format(
                                                    pdb_stem.mids[0], target_stem.mids[0], d_start)
    assert d_start<0.1, "{Distance between stem ends {} and {} is too big: {}".format(
                                                    pdb_stem.mids[1], target_stem.mids[1], d_end)
    tw1_polar_pdb    = ftuv.spherical_cartesian_to_polar(pdb_stem.twists[0])
    tw1_polar_target = ftuv.spherical_cartesian_to_polar(target_stem.twists[0])
    d_twist_u = abs(tw1_polar_pdb[1]-tw1_polar_target[1])
    d_twist_v = abs(tw1_polar_pdb[2]-tw1_polar_target[2])
    assert d_twist_u<0.01, "Deviation of twist angle u too big: {:f}".format(d_twist_u)
    assert d_twist_v<0.01, "Deviation of twist angle v too big: {:f}".format(d_twist_u)
    return True

def translate_chain(chains, translation):
    '''
    Translate all of the atoms in a chain by a certain amount.

    :param chains: A dict chain_id:Bio.PDB.Chain instance to be translated.
    :translation: A vector indicating the direction of the translation.
    '''

    for chain in chains.values():
        atoms = bpdb.Selection.unfold_entities(chain, 'A')
        for atom in atoms:
            atom.transform(ftuv.identity_matrix, translation)

def rotate_chain(chains, rot_mat, offset):
    '''
    Move according to rot_mat for the position of offset.

    :param chains: A dict chain_id:Bio.PDB.Chain instance.
    :param rot_mat: A left_multiplying rotation_matrix.
    :param offset: The position from which to do the rotation.
    '''
    new_coords = []
    for chain in chains.values():

        atoms = bpdb.Selection.unfold_entities(chain, 'A')

        for atom in atoms:
            #atom.transform(ftuv.identity_matrix, -offset)
            assert ftuv.magnitude(atom.coord - offset) < ftuv.magnitude(atom.coord)
            atom.coord -= offset
            new_coords.append(atom.coord)
            atom.transform(rot_mat, offset)
        dev_from_cent = ftuv.magnitude(np.sum(new_coords, axis=0)/len(new_coords))
        if dev_from_cent>1:
            log.warning("{} not close to zero".format(dev_from_cent))



def get_stem_rotation_matrix(stem, stem2, use_average_method=False):
    """
    :param stem: The first StemModel
    :param stem2: The second StemModel

    :retuirns: A RotationMatrix.
               Use stem1.vec()*rotMat to rotate stem1 onto stem2
               Use rotMat*stem2.vec() to rotate stem2 onto stem1
    """
    #twist1 = (stem.twists[0] + stem.twists[1]) / 2.

    if not use_average_method:
        twist1 = stem.twists[0]
        twist2 = stem2.twists[0]

    else:
        twist1 = ftug.virtual_res_3d_pos_core(stem.mids, stem.twists, 2, 4)[1]
        twist2 = ftug.virtual_res_3d_pos_core(stem2.mids, stem2.twists, 2, 4)[1]

    return ftuv.get_double_alignment_matrix((stem.vec(), twist1),(stem2.vec(), twist2))
    # get normalvector to stem and twist.
    comp1 = np.cross(stem.vec(), twist1)

    # rotate around the first stem by t degrees
    rot_mat1 = ftuv.rotation_matrix(stem.vec(), t)
    rot_mat2 = ftuv.rotation_matrix(twist1, u - math.pi/2)
    rot_mat3 = ftuv.rotation_matrix(comp1, v)

    rot_mat4 = np.dot(rot_mat3, np.dot(rot_mat2, rot_mat1))

    return rot_mat4

def insert_element(cg_to, cg_from, elem_to, elem_from,
                   chains_to, chains_from):
    '''
    Take an element (elem_from) from one dict of chains (chains_from, cg_from) and
    insert it on the new chain while aligning on the adjoining elements.

    The neighboring elements need to be present in chain_to in order
    for the next element to be aligned to their starting and ending
    positions.

    The dimensions and type of elem_to and elem_from need to be identical.

    This method aligns the flanking base pairs on both ends
    (except for 3' and 5' elements) of the fragment with the respective
    base-pairs in the stem-scaffold. This means that there will be
    equally big breaks in the chain on both sides of the fragment.


    @param cg_to: The coarse-grain representation of the target chain
    @param cg_from: The coarse-grain representation of the source chain
    @param elem_to: The element to replace
    @param elem_from: The source element
    @param chains_to: A dict chainid:chain. The chains to graft onto
    @param chains_from: A dict chainid:chain. The chains to excise from
    '''

    assert elem_from[0]==elem_to[0]
    # The define of the loop with adjacent nucleotides (if present) in both cgs
    define_a_to = cg_to.define_a(elem_to)
    define_a_from = cg_from.define_a(elem_from)
    assert len(define_a_to) == len(define_a_from)
    # The defines translated to seq_ids.
    closing_bps_to = []
    closing_bps_from = []
    for nt in define_a_to:
        closing_bps_to.append(cg_to.seq_ids[nt-1])
    for nt in define_a_from:
        closing_bps_from.append(cg_from.seq_ids[nt-1])
    # Seq_ids of all nucleotides in the loop that will be inserted
    seq_ids_a_from = []
    for i in range(0, len(define_a_from), 2):
        for nt in range(define_a_from[i], define_a_from[i+1]+1):
            seq_ids_a_from.append(cg_from.seq_ids[nt-1])
    #The loop fragment to insert in a dict {chain_id:chain}
    try:
        pdb_fragment_to_insert = ftup.extract_subchains_from_seq_ids(chains_from, seq_ids_a_from)
    except:
        log.error("Could not extract fragment %s from pdb: "
                  " At least one of the seq_ids %s not found."
                  " Chains are %s", elem_from, seq_ids_a_from, chains_from.keys())
        raise

    # A list of tuples (seq_id_from, seq_id_to) for the nucleotides
    # that will be used for alignment.
    log.info("Closing_bps _from are %s", closing_bps_from)
    alignment_positions = []
    if elem_from[0]=="t": #Use only left part of define
        alignment_positions.append((closing_bps_from[0], closing_bps_to[0]))
    elif elem_from[0]=="f": #Use only right part of define
        alignment_positions.append((closing_bps_from[1], closing_bps_to[1]))
    else: #Use all flanking nucleotides
        assert elem_from[0]!="s", "No stems allowed in insert_element"
        for i in range(len(closing_bps_from)):
            alignment_positions.append((closing_bps_from[i], closing_bps_to[i]))
    align_on_nucleotides(chains_from, chains_to, alignment_positions)

    #The defines and seq_ids WITHOUT adjacent elements
    define_to = cg_to.defines[elem_to]
    define_from = cg_from.defines[elem_from]
    seq_ids_to = []
    seq_ids_from = []
    for i in range(0, len(define_from), 2):
        for nt in range(define_from[i], define_from[i+1]+1):
            seq_ids_from.append(cg_from.seq_ids[nt-1])
        for nt in range(define_to[i], define_to[i+1]+1):
            seq_ids_to.append(cg_to.seq_ids[nt-1])
    assert len(seq_ids_to)==len(seq_ids_from)
    # Now copy the residues from the fragment chain to the scaffold chain.
    for i in range(len(seq_ids_from)):
        resid_from = seq_ids_from[i]
        resid_to   = seq_ids_to[i]

        residue = chains_from[resid_from.chain][resid_from.resid]
        #Change the resid to the target
        residue.id = resid_to.resid
        if resid_to.chain not in chains_to:
            log.info("Adding chain with id %r for residue %r", resid_to.chain, resid_to)
            chains_to[resid_to.chain] =  bpdb.Chain.Chain(resid_to.chain)
        #Now, add the residue to the target chain
        chains_to[resid_to.chain].add(residue)

def align_on_nucleotides(chains_fragment, chains_scaffold, nucleotides):
    """
    Translate and rotate pdb-fragment to optimally match
    the scaffold at the specified nucleotides.

    :param chains_scaffold: A dict {chain_id:chain} that contains the stems
                            where the fragment should be inserted.
    :param chains_fragment: A dict {chain_id:chain} that contains the fragment
                            to be inserted (with adjacent nucleotides)
    :param nucleotides: A list of tuples (seq_id_fragment, seq_id_scaffold)
                        The two seq_ids in each tuple spould refer to the same
                        nucleotide (the fragments adjacent nts), once in the
                        fragment, once in the scaffold.
    """
    # The point-clouds that will be aligned.
    points_fragment = []
    points_scaffold = []
    log.info("nts %s", nucleotides)
    for res_frag, res_scaf in nucleotides:
        #log.debug("res_frag %s res_scaf %s", res_frag, res_scaf)
        residue_frag = chains_fragment[res_frag.chain][res_frag.resid]
        residue_scaf = chains_scaffold[res_scaf.chain][res_scaf.resid]
        for atom_label in ["C4'", "C3'", "O3'"]:
            points_fragment.append(residue_frag[atom_label].coord)
            points_scaffold.append(residue_scaf[atom_label].coord)
    points_fragment = np.asarray(points_fragment)
    points_scaffold = np.asarray(points_scaffold)

    centroid_frag = ftuv.get_vector_centroid(points_fragment)
    centroid_scaf = ftuv.get_vector_centroid(points_scaffold)

    points_fragment -= centroid_frag
    points_scaffold -= centroid_scaf

    sup = brmsd.optimal_superposition(points_fragment, points_scaffold)

    for chain in chains_fragment.values():
        for atom in chain.get_atoms():
            atom.transform(np.eye(3,3), -centroid_frag)
            atom.transform(sup, centroid_scaf)

# used for Validation
def _compare_cg_chains_partial(cg, chains):
    """
    :param cg: The CoarseGrainRNA
    :param chains: A dict {chain_id: Chain}
    """
    for chain in chains.values():
        for res in chain.get_residues():
            resid = fgb.RESID(chain.id, res.id)
            pdb_coords = res["C1'"].coord
            cg_coords = cg.virtual_atoms(cg.seq_ids.index(resid)+1)["C1'"]
            if ftuv.magnitude(pdb_coords-cg_coords)>2:
                log.error("Residue %s, C1' coords %s do not "
                          "match the cg-coords (virtual atom) %s by %f", resid, pdb_coords,
                          cg_coords, ftuv.magnitude(pdb_coords-cg_coords))


def align_residues(res_dir, res_ref):
    '''
    Orient res_ref so that it points in the same direction
    as res_dir.

    :param res_dir: The residue indicating the direction
    :param res_ref: The reference residue to be rotated
    :return res: A residue with the atoms of res_ref pointing in the direction of res_dir
    '''
    #TODO: BT: In the future we might align based on ALL non-sidechain atoms using optimal_superposition.
    #     If we stick to 3 reference atoms, we could use ftuv.get_double_alignment_matrix instead.

    av = { 'U': ['N1', "C1'", "C2'"], 'C': ['N1', "C1'", "C2'"], 'A': ['N9', "C1'", "C2'"], 'G': ['N9', "C1'", "C2'"],
           'rU': ['N1', "C1'", "C2'"], 'rC': ['N1', "C1'", "C2'"], 'rA': ['N9', "C1'", "C2'"], 'rG': ['N9', "C1'", "C2'"] }

    dv = av[res_dir.resname.strip()]
    rv = av[res_ref.resname.strip()]

    dir_points = np.array([res_dir[v].get_vector().get_array() for v in dv])
    ref_points = np.array([res_ref[v].get_vector().get_array() for v in rv])

    dir_centroid = cuv.get_vector_centroid(dir_points)
    ref_centroid = cuv.get_vector_centroid(ref_points)

    sup = brmsd.optimal_superposition(ref_points - ref_centroid, dir_points - dir_centroid)
    new_res = res_ref.copy()

    for atom in new_res:
        atom.transform(np.eye(3,3), -ref_centroid)
        atom.transform(sup, dir_centroid)

    return new_res

def replace_bases(chains, cg):
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

    if fbc.Configuration.template_residues is None:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            fbc.Configuration.template_residues = bpdb.PDBParser().get_structure('t', conf.Configuration.template_residue_fn)
    s1 = fbc.Configuration.template_residues
    tchain = list(s1.get_chains())[0]


    templates = { 'A': tchain[1], 'C': tchain[2], 'G': tchain[3], 'U': tchain[4]}

    for chain_name, chain in chains.items():
        residues = chain.get_list()

        for residue in residues:
            #num = ress[i].id[1]
            old_name = residue.resname.strip()
            cg_seq_num = cg.seq_ids.index(fgb.RESID(chain = chain_name, resid = residue.id))
            target_name = cg.seq[cg_seq_num]
            if target_name == old_name:
                #Don't replace the correct residues
                log.debug("Not replacing %s by %s", old_name, target_name )
                continue
            log.debug("Replacing %s by %s", old_name, target_name)
            ref_res = templates[target_name]
            new_res = align_residues(residue, ref_res)
            sca = ftup.side_chain_atoms[old_name]
            for aname in sca:
                residue.detach_child(aname)

            sca = ftup.side_chain_atoms[target_name]
            for aname in sca:
                residue.add(new_res[aname])

            residue.resname = target_name

def reorder_residues(chains, cg):
    '''
    Reorder the nucleotides in the chain's list so that they match the order
    in the cg representation.

    :param chains: A dict {chain_id:Bio.PDB.Chain}
    :param cg: A coarse grain representation
    '''
    for chain_name, chain in chains.items():
        chain.child_list.sort(key=lambda x: cg.seq_ids.index(fgb.RESID(chain = chain_name, resid = x.id)))
