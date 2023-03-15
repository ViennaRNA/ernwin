'''
Convert a bulge graph to the pdb-format.

'''


from __future__ import absolute_import
import multiprocessing as mp
import warnings
import glob
from collections import defaultdict

import itertools as it
import fess.builder.models as models
import fess.builder.config as fbc
from fess import data_file
import os.path as op

from forgi.graph.residue import RESID
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.pdb as ftup
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.vector as cuv
from six.moves import map
from six.moves import range
ftuv = cuv
import forgi.utilities.debug as fud
import forgi.utilities.stuff as fus

import forgi.threedee.model.similarity as brmsd
import forgi.threedee.model.stats as ftms

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

from logging_exceptions import log_to_exception
import logging
log=logging.getLogger(__name__)

import numpy as np


class Reconstructor(object):

    def __init__(self, pdb_library_path, cg_library_path, server=False, library_directory=None):
        """
        :param server: If True, assume a lot of memory is available
        """
        self.cg_library_path = cg_library_path
        self.pdb_library_path = pdb_library_path
        self.stem_use_average_method = True
        self._loaded_cgs = {}
        self._loaded_pdbs = {}
        self.store=server
        self.library_directory=library_directory
    def reconstruct(self, sm):
        '''
        Re-construct a full-atom model from a coarse-grain model.

        :param sm: The Spatial Model to reconstruct
        :returns: A dictionary {chain_id : Bio.PDB.Chain.Chain}
        '''
        # First, reconstruct the stems, which will be used as scaffold for the bulges
        chains = {}
        num_stems = len(list(sm.bg.stem_iterator()))
        for i, stem in enumerate(sm.bg.stem_iterator()):
            self._reconstruct_stem(sm, stem, chains)
            log.info("{} out of {} stems reconstructed".format(i+1, num_stems))
        #Some validation
        try:
            _compare_cg_chains_partial(sm.bg, chains)
        except Exception:
            log.exception("Validation could not be performed due to the following error!")


        # Some more verification
        for res in sm.bg.seq._seqids:
            if sm.bg.get_elem(res)[0]=="s":
                try:
                    chains[res.chain][res.resid]
                except LookupError as e:
                    log.error("STEM Residue missing in reconstructed PDB: %s (KeyError %s)", res, e)
                    log.error("It belongs to %s", sm.bg.get_elem(res))


        #ftup.output_multiple_chains(chains.values(), "reconstr_stems.pdb")
        '''# Some more validation. Useless unless searching for a bug
        chains_in_file = ftup.get_all_chains("reconstr_stems.pdb")
        chains_in_file = { c.id:c for c in chains_in_file }
        for c in chains:
            r = ftup.pdb_rmsd(chains[c], chains_in_file[c], superimpose = False)[1]
            assert r<0.001, "r={} for chain {}".format(r, c)
        '''
        gaps_to_mend = []
        for loop in sm.bg.defines:
            if loop[0]=="s": continue
            gaps_to_mend.extend(self._reconstruct_with_fragment(chains, sm, loop))
        #ftup.output_multiple_chains(chains.values(), "reconstr_all.pdb")

        replace_bases(chains, sm.bg)
        #ftup.output_multiple_chains(chains.values(), "base_replaced.pdb")
        reorder_residues(chains, sm.bg)
        chains = mend_breakpoints(chains, gaps_to_mend)

        #ftup.output_multiple_chains(chains.values(), "mended.pdb")
        return chains

    def _get_fragment(self, stat, sm):
        key = stat.pdb_name+"__def_"+"-".join(map(str,stat.define))
        new_fragment=False
        try:
            fragment,_,_ = ftup.get_all_chains(op.join(self.library_directory, key[2:4], key+".cif"),
                                             no_annotation=True)
        except Exception:
            cg, chains = self._get_source_cg_and_chain(stat, sm)
            new_fragment=True
        else:
            fragment = {c.id:c for c in fragment}
            log.debug("Used stored fragment for %s", key)
            pdb_basename = stat.pdb_name.split(":")[0]
            cg_filename = op.expanduser(op.join(self.cg_library_path, pdb_basename+".cg"))
            cg = self.get_cg(cg_filename)  #The cg with the template
        try:
            elem = cg.get_node_from_residue_num(stat.define[0])
        except Exception:
            log.error("stat %s with define %s", stat, stat.define)
            raise
        if stat.define != cg.defines[elem]:
            err = ValueError("The CG files where the stats where extracted and "
                             "the cg file used for reconstruction are not consistent!")
            with log_to_exception(log, err):
                log.error("%s != %s for element %s (%s)", stat.define, cg.defines[elem], elem, stat.pdb_name)
            raise err
        if new_fragment:
            fragment  = ftup.extract_subchains_from_seq_ids(chains,
                            cg.define_residue_num_iterator(elem, seq_ids=True,
                                                           adjacent=(elem[0]!="s")))
            if self.library_directory is not None:
                log.debug("Storing newly-created fragment for %s", key)
                import distutils.dir_util
                distutils.dir_util.mkpath(op.join(self.library_directory, key[2:4]))
                ftup.output_multiple_chains(list(fragment.values()),
                                            op.join(self.library_directory, key[2:4], key+".cif"), "cif")
        return cg, elem, fragment

    def _get_source_cg_and_chain(self, stat, sm):
        """
        Load the fragment defined in the stat from the fragment library as pdb and cg.

        :param stat: The forgi.threedee.model.stats.StemStat or ftms.AngleStat or ftms.LoopStat object.
        :param sm: The SpatialModel to reconstruct. Used, if it contains stats not sampled but loaded directly.
        """
        stat_name = stat.pdb_name
        if stat_name == sm.bg.name and sm.bg.chains:
            return sm.bg, sm.bg.chains

        pdb_basename = stat_name.split(":")[0]
        pdb_filename = op.expanduser(op.join(self.pdb_library_path, "_".join(pdb_basename.split("_")[:-1])+".pdb"))
        cg_filename = op.expanduser(op.join(self.cg_library_path, pdb_basename+".cg"))
        #Make sure the files exist.
        try:
            try:
                with open(pdb_filename): pass
            except IOError:
                pdb_filename = pdb_filename.rstrip(".pdb")+".cif"
                with open(pdb_filename): pass
            with open(cg_filename): pass
        except Exception as e:
            with log_to_exception(log, e):
                log.error("Failed to open files for stat %s", stat.pdb_name)
            raise
        log.debug("Opening cg-file %s to extract stat %s", cg_filename, stat.pdb_name)
        cg = self.get_cg(cg_filename)  #The cg with the template
        chains = self.get_pdb(pdb_filename, store = self.store)

        return cg, chains

    def get_cg(self, fn):
        if fn not in self._loaded_cgs:
            cg = ftmc.CoarseGrainRNA.from_bg_file(fn)
            return cg
        return self._loaded_cgs[fn]

    def get_pdb(self, fn, store=False):
        if fn not in self._loaded_pdbs:
            chains, missing_residues, _ = ftup.get_all_chains(fn, no_annotation=True)
            new_chains = []
            for chain in chains:
                chain, modifications = ftup.clean_chain(chain)
                new_chains.append(chain)
            chains = {c.id:c for c in new_chains}
            if store:
                self._loaded_pdbs[fn]=chains
            return chains
        return self._loaded_pdbs[fn]

    def _reconstruct_stem(self, sm, stem_name, new_chains):
        '''
        Reconstruct a particular stem.

        :param sm: The SpatialModel that should be translated to PDB
        :param stem_name: The name of the stem to be reconstructed.
        :param new_chains: A dict chainid:chain that will be filled with the reconstructed model
        '''
        assert stem_name[0]=="s"
        stem = sm.stems[stem_name]
        stem_stat = sm.elem_defs[stem_name]
        orig_def = sm.bg.defines[stem_name]
        cg_orig = sm.bg
        cg, sd, chains = self._get_fragment(stem_stat, sm)

        _align_chain_to_stem(cg, chains, sd, stem, self.stem_use_average_method)

        for i in range(stem_stat.bp_length):
            for strand in range(2):
                target_resid = cg_orig.seq.to_resid(orig_def[strand*2] + i)
                source_resid = cg.seq.to_resid(stem_stat.define[strand*2] + i)
                residue = chains[source_resid.chain][source_resid.resid]
                #Change the resid to the target
                residue.parent=None
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
            s1,s2 = sorted(sm.bg.edges[ld], key=lambda x: sm.bg.defines[x][0])
            log.warning("%s Not part of MST. Trying to extract stat from SM", ld)
            angle_stat=sm.bg.get_bulge_angle_stats_core(ld,(s1,s2))
            log.debug("angle stat is %s", angle_stat)
        if not angle_stat.define:
            # A zero-length multiloop. We do not need to insert any fragment.
            return []
        cg_from, elem_from, chains_from = self._get_fragment(angle_stat, sm)
        cg_to = sm.bg
        chains_to = chains
        elem_to = ld

        if isinstance(angle_stat, ftms.AngleStat):
            angle_type = angle_stat.ang_type
        else:            # LoopStat has no angle type
            angle_type = 0
        return insert_element(cg_to, cg_from, elem_to, elem_from,
                              chains_to, chains_from, angle_type)


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
    except Exception as e:
        with log_to_exception(log, e):
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
    if d_twist_u>0.01:
        log.info("Deviation of twist angle u too big for %s: %s", elem_name, d_twist_u)
    if d_twist_v>0.01:
        log.info("Deviation of twist angle v too big for %s: %s", elem_name, d_twist_v)
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
            #assert ftuv.magnitude(atom.coord - offset) < ftuv.magnitude(atom.coord), "{}!<{}".format(ftuv.magnitude(atom.coord - offset), ftuv.magnitude(atom.coord))
            atom.coord -= offset
            new_coords.append(atom.coord)
            atom.transform(rot_mat, offset)
        dev_from_cent = ftuv.magnitude(np.sum(new_coords, axis=0)/len(new_coords))
        if dev_from_cent>5:
            log.info("{} not close to zero".format(dev_from_cent))



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

def iter_resids_between(chain, from_, to, include_start, include_end):
    log.debug("Iter resids vewtween %s and %s from %s", from_, to, [x.id for x in chain.child_list])
    found = False
    for res in chain:
        if res.id == from_:
            found=True
            if include_start:
                log.debug("yielding start %s", res.id)
                yield RESID(chain.id, res.id)
        elif res.id == to:
            if include_end:
                log.debug("yielding end %s", res.id)
                yield RESID(chain.id, res.id)
            return
        elif found:
            log.debug("yielding %s", res.id)
            yield RESID(chain.id, res.id)


def insert_element(cg_to, cg_from, elem_to, elem_from,
                   chains_to, chains_from, angle_type):
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


    :param cg_to: The coarse-grain representation of the target chain
    :param cg_from: The coarse-grain representation of the source chain
    :param elem_to: The element to replace
    :param elem_from: The source element
    :param chains_to: A dict chainid:chain. The chains to graft onto
    :param chains_from: A dict chainid:chain. The chains to excise from

    :returns: a list of tuples containing gaps to mend
    '''
    log.info("Inserting element %s",elem_to)
    assert elem_from[0]==elem_to[0], "{}[0]!={}[0]".format(elem_from, elem_to)
    # The define of the loop with adjacent nucleotides (if present) in both cgs
    define_a_to = cg_to.define_a(elem_to)
    define_a_from = cg_from.define_a(elem_from)
    assert len(define_a_to) == len(define_a_from)
    nt_in_define_from=[x in cg_from.defines[elem_from] for x in define_a_from]

    # The defines translated to seq_ids.
    closing_bps_to = []
    closing_bps_from = []

    log.debug("Angle type is %s", angle_type)
    for nt in define_a_to:
        closing_bps_to.append(cg_to.seq.to_resid(nt))
    for nt in define_a_from:
        closing_bps_from.append(cg_from.seq.to_resid(nt))
    # Seq_ids of all nucleotides in the loop that will be inserted
    seq_ids_a_from = []
    for i in range(0, len(define_a_from), 2):
        for nt in range(define_a_from[i], define_a_from[i+1]+1):
            seq_ids_a_from.append(cg_from.seq.to_resid(nt))
    log.debug("seqids_a from %s", seq_ids_a_from)
    #The loop fragment to insert in a dict {chain_id:chain}
    try:
        chains_from = ftup.extract_subchains_from_seq_ids(chains_from, seq_ids_a_from)
    except Exception as e:
        with log_to_exception(log, e):
            log.error("Could not extract fragment %s from pdb: "
                      " At least one of the seq_ids %s not found."
                      " Chains are %s", elem_from, seq_ids_a_from, list(chains_from.keys()))
        raise

    # A list of tuples (seq_id_from, seq_id_to) for the nucleotides
    # that will be used for alignment.
    log.debug("Closing_bps _from are %s", closing_bps_from)
    alignment_positions = []
    assert elem_from[0]!="s", "No stems allowed in insert_element"
    if elem_from[0]=="f":
        alignment_positions.append((closing_bps_from[1], closing_bps_to[1]))
    elif elem_from[0]=="t":
        alignment_positions.append((closing_bps_from[0], closing_bps_to[0]))
    else:
        for i in range(len(closing_bps_from)): #
            alignment_positions.append((closing_bps_from[i], closing_bps_to[i]))

    log.debug("Calling align_on_nucleotides for %s", elem_to)
    align_on_nucleotides(chains_from, chains_to, alignment_positions)

    #The defines and seq_ids WITHOUT adjacent elements
    define_to = cg_to.defines[elem_to]
    define_from = cg_from.defines[elem_from]
    no_moderna = False
    if len(define_from)!=len(define_to):
        log.warning("Inconsistent defines: {} and {} for {}. Using ModeRNA fragment instead.".format(define_from, define_to, elem_to))
        target_seqs = cg_to.get_define_seq_str(elem_to) # One or two strands
        for i, target_seq in enumerate(target_seqs):
            if closing_bps_from[2*i].chain!=closing_bps_from[2*i+1].chain:
                raise NotImplementedError("TODO")
            try:
                mod_chain = use_moderna_fragment(chains_from[closing_bps_from[2*i].chain], target_seq,
                                             closing_bps_from[2*i], closing_bps_from[2*i+1] )
            except:
                no_moderna=True
            else:
                chains_from[seq_ids_a_from[0].chain]= mod_chain
    elif cg_to.element_length(elem_to) != cg_from.element_length(elem_from):
        log.warning("%s not consistent with %s: Missing residues for %s", define_from, define_to, elem_to)
        log.warning("%s has different len than %s for angle type %s", define_from, define_to, angle_type)
        if define_to[1]-define_to[0]>define_from[1]-define_from[0]:
            # Apply an indel on the left side
            if closing_bps_from[0].chain!=closing_bps_from[1].chain:
                raise NotImplementedError("TODO")
            target_seq = cg_to.get_define_seq_str(elem_to)[0] # Forward strand
            try:
                mod_chain = use_moderna_fragment(chains_from[closing_bps_from[0].chain], target_seq,
                                                 closing_bps_from[0], closing_bps_from[1])
            except:
                no_moderna = True
            else:
                chains_from[seq_ids_a_from[0].chain]= mod_chain
        else:
            log.warning("Reconstruction for %s will contain nucleotides missing from the input structure", elem_to)
            if elem_to[0]=="h" and define_to[1]-define_to[0]<3:
                raise NotImplementedError("Reconstruction of hairpins with length ,3 is not supportet!")
            #target_seq = cg_to.get_define_seq_str(elem_to)[0] # Forward strand
            #mod_chain = use_moderna_fragment(chains_from[closing_bps_from[0].chain], target_seq,
            #                                     closing_bps_from[0], closing_bps_from[1])
            #chains_from[seq_ids_a_from[0].chain]= mod_chain
    seq_ids_to = []
    for i in range(0, len(define_to), 2):
        seq_ids_to.append([])
        for seq_id in cg_to.seq.with_missing.iter_resids(cg_to.seq.to_resid(define_to[i]),cg_to.seq.to_resid(define_to[i+1])):
            seq_ids_to[-1].append(seq_id)
    seq_ids_from = []
    # Now append first strand to seq_ids_from
    assert closing_bps_from[0].chain == closing_bps_from[1].chain
    log.debug("define_from=%s, define_to=%s, nt_in_define=%s", define_from, define_to, nt_in_define_from)
    if True: #closing_bps_from[0].resid<closing_bps_from[1].resid:
        s = list(iter_resids_between(chains_from[closing_bps_from[0].chain], closing_bps_from[0].resid, closing_bps_from[1].resid, nt_in_define_from[0], nt_in_define_from[1]))
    else:
        s = list(iter_resids_between(chains_from[closing_bps_from[0].chain], closing_bps_from[1].resid, closing_bps_from[0].resid, nt_in_define_from[1], nt_in_define_from[0]))
        s[0].reverse()
    if s:
        seq_ids_from.append(s)

    if len(closing_bps_from)>2:
        assert closing_bps_from[2].chain == closing_bps_from[3].chain
        if closing_bps_from[2].resid<closing_bps_from[3].resid:
            s = (list(iter_resids_between(chains_from[closing_bps_from[2].chain], closing_bps_from[2].resid, closing_bps_from[3].resid, nt_in_define_from[2], nt_in_define_from[3])))
        else:
            s = list(iter_resids_between(chains_from[closing_bps_from[2].chain], closing_bps_from[3].resid, closing_bps_from[2].resid, nt_in_define_from[3], nt_in_define_from[2]))
            s.reverse()
        if s:
            seq_ids_from.append(s)



    log.info("Fragment %s", seq_ids_from)
    log.info("Target %s", seq_ids_to)
    #if not no_moderna:
    #    assert len(seq_ids_from[0]) == len(seq_ids_to[0]), "Unequal length for {}: {} {}".format(elem_to, seq_ids_from, seq_ids_to)
    #    if len(seq_ids_to)>1:
    #        assert len(seq_ids_from[1]) == len(seq_ids_to[1]), "Unequal length for {}: {} {}".format(elem_to, seq_ids_from, seq_ids_to)
    log.debug("Copying %s to %s for %s", seq_ids_from, seq_ids_to, elem_to)
    # Now copy the residues from the fragment chain to the scaffold chain.
    lastres=[None, None]
    for a in range(len(seq_ids_to)):
        for i in range(len(seq_ids_to[a])):
            try:
                resid_from = seq_ids_from[a][i]
            except IndexError:
                lastres[a]= seq_ids_to[a][i-1]
                break
            resid_to   = seq_ids_to[a][i]
            residue = chains_from[resid_from.chain][resid_from.resid]
            #Change the resid to the target
            residue.parent = None
            residue.id = resid_to.resid
            if resid_to.chain not in chains_to:
                log.info("Adding chain with id %r for residue %r", resid_to.chain, resid_to)
                chains_to[resid_to.chain] =  bpdb.Chain.Chain(resid_to.chain)
            #Now, add the residue to the target chain
            chains_to[resid_to.chain].add(residue)
    # Now we need to mend gaps created by imperfect alignment.
    gaps_to_mend = []
    if elem_from[0]!="f":
        log.debug("To mend: %s %s ", cg_to.seq.to_resid(define_a_to[0]), cg_to.seq.to_resid(define_a_to[0]+1))
        gaps_to_mend.append( [ cg_to.seq.to_resid(define_a_to[0]),
                               cg_to.seq.to_resid(define_a_to[0]+1) ] )
        d = gap_length(chains_to, cg_to.seq.to_resid(define_a_to[0]),  cg_to.seq.to_resid(define_a_to[0]+1))
        log.debug("Elem {}: dist {} - {} is {}".format(elem_to, define_a_to[0], define_a_to[0]+1, d))
    if elem_from[0]!="t":
        if lastres[0] is not None:
            r=lastres[0]
        else:
            r=cg_to.seq.to_resid(define_a_to[1]-1)
        gaps_to_mend.append( [ r, cg_to.seq.to_resid(define_a_to[1]) ] )
        d = gap_length(chains_to, r,  cg_to.seq.to_resid(define_a_to[1]))
        log.debug("Elem {}: dist {} - {} is {}".format(elem_to, define_a_to[1], define_a_to[1]-1, d))

    if elem_from[0]=="i":
        gaps_to_mend.append( [ cg_to.seq.to_resid(define_a_to[2]),
                               cg_to.seq.to_resid(define_a_to[2]+1) ] )
        if lastres[1] is not None:
            r=lastres[1]
        else:
            r=cg_to.seq.to_resid(define_a_to[3]-1)
        gaps_to_mend.append( [ r,
                               cg_to.seq.to_resid(define_a_to[3]) ] )
    log.debug("To mend %s", gaps_to_mend)
    return gaps_to_mend

def use_moderna_fragment(chain, target_seq, resid1, resid2):
    try:
        import moderna
    except ImportError:
        raise RuntimeError("ModeRNA is required for processing indels!")
    mod_model = moderna.load_model(chain, data_type="chain")
    res5p = resid_to_moderna(resid1)
    res3p = resid_to_moderna(resid2)
    moderna.apply_indel(mod_model, res5p,
                        res3p,
                        str(target_seq)
                        )
    # Back to PDB
    pdb_model = mod_model.get_structure()[0]
    mod_chain, = pdb_model.child_list
    return mod_chain


def gap_length(chains, resid1, resid2):
    try:
        res1 = chains[resid1.chain][resid1.resid]
        res2 = chains[resid2.chain][resid2.resid]
    except KeyError as e:
        log.exception("Cannot calculate gap length")
        return None
    try:
        d = ftuv.vec_distance(res1["P"].coord, res2["P"].coord)
    except KeyError as e:
        log.exception("%s not in %s or %s",e,res1.child_dict, res2.child_dict)
        return None
    return d

def resid_to_moderna(resid):
    """
    :param resid: fgr.RESID
    :returns: A string of the form "123" or "123B"
    """
    res = resid.resid
    res = (str(res[1])+res[2]).strip()
    return res

def mend_breakpoints(chains, gap):
    """
    :param gap: A list of res_ids, which can be moved to mend the gap.
    """
    #raise NotImplementedError("Error")
    try:
        import moderna
    except ImportError:
        warnings.warn("Cannot mend gaps in sequence, because ModeRNA is not installed!")
        return chains
    mod_models = {}
    with fus.make_temp_directory() as tmpdir:
        log.info("Writing chains %s", list(chains.values()))

        #ftup.output_multiple_chains(chains.values(), op.join(tmpdir, "tmp.pdb"))
        for g in gap:
            if g[0].chain != g[1].chain:
                log.warning("Not mending gap between multiple chains: %s and %s",
                            g[0], g[1])
                continue
            if g[0].chain not in mod_models:
                try:
                    mod_models[g[0].chain] =  moderna.load_model(chains[g[0].chain], data_type="chain")#moderna.load_model(op.join(tmpdir, "tmp.pdb"), g[0].chain)
                except Exception as e:
                    with log_to_exception(log, e):
                        log.error("g is %s, g[0] is %s, g[0].chain is %s", g, g[0], g[0].chain)
                        log.error("chains is %s", chains)
                    raise
            moderna.fix_backbone(mod_models[g[0].chain], resid_to_moderna(g[0]), resid_to_moderna(g[1]))
            #moderna.write_model(mod_models[g[0].chain], op.join(tmpdir, "tmp.pdb"))
        #for chain_id, model in mod_models.items():
        #    moderna.write_model(model,  op.join(tmpdir, "mended_{}.pdb".format(chain_id)))
        #Load back to Biopython
        mended_chains = {}
        for chain_id in chains.keys():
            if chain_id in mod_models:
                mended_chains[chain_id] = mod_models[chain_id] #Mod models are chain subclasses anyway
                log.info("Mended:", mended_chains)
                mended_chains[chain_id].id = chain_id
            else:
                mended_chains[chain_id] = chains[chain_id]
    log.info("mended_chains: %s", mended_chains)
    # Moderna may replace modified residues with "UNK" for unknown or otherrwise change the code.
    # We have to replace them back.
    for chain_id in chains:
        for res in mended_chains[chain_id]:
            changed = False
            for o_res in chains[chain_id]:
                if o_res.id[1:]==res.id[1:]:
                    log.debug("Changing Moderna residue %s to %s", res, o_res)
                    assert not changed #Only one residue per number+icode
                    res.id = o_res.id
                    res.resname = o_res.resname
                    log.debug("Moderna residue now %s", res)
                    changed = True
    # Convert back from ModeRNA to Biopython
    out_chains={}
    for k,v in mended_chains.items():
        s = v.get_structure()[0]
        log.error("%s, %s %s", k, s, s.child_dict)
        assert len(s.child_list)==1
        out_chains[k]=s.child_list[0]
        out_chains[k].id=k
    return out_chains
    #{k:v.get_structure()[0][k] for k,v in mended_chains.items()}

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
        found=0
        for atom_label in ["C5'", "C4'", "C3'", "O3'"]:
            try:
                frag_coords = residue_frag[atom_label].coord
                scaf_coords = residue_scaf[atom_label].coord
            except KeyError:
                pass
            else:
                found+=1
                points_fragment.append(residue_frag[atom_label].coord)
                points_scaffold.append(residue_scaf[atom_label].coord)
        if found<3:
            raise ValueError("Cannot align_on_nucleotides: Too many missing atoms")
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
            log.debug("Resid %s", resid)
            virt_a = cg.virtual_atoms(cg.seq.to_integer(resid))
            try:
                cg_coords = virt_a["C1'"]
            except KeyError:
                log.error("virtual atom coordinates for %s only for %s", resid, list(virt_a.keys()))
                raise
            if ftuv.magnitude(pdb_coords-cg_coords)>4:
                log.warning("Residue %s, C1' coords %s do not "
                          "match the cg-coords (virtual atom) %s by %f", resid, pdb_coords,
                          cg_coords, ftuv.magnitude(pdb_coords-cg_coords))


def align_residues(res_dir, res_ref, on=None):
    '''
    Orient res_ref so that it points in the same direction
    as res_dir.

    :param res_dir: The residue indicating the direction
    :param res_ref: The reference residue to be rotated
    :return res: A residue with the atoms of res_ref pointing in the direction of res_dir
    '''
    #TODO: BT: In the future we might align based on ALL non-sidechain atoms
    #     using optimal_superposition.
    #     If we stick to 3 reference atoms, we could use
    #     ftuv.get_double_alignment_matrix instead.

    if on is None:
        av = { 'U': ['N1', "C1'", "C2'"], 'C': ['N1', "C1'", "C2'"], 'A': ['N9', "C1'", "C2'"], 'G': ['N9', "C1'", "C2'"],
              'rU': ['N1', "C1'", "C2'"], 'rC': ['N1', "C1'", "C2'"], 'rA': ['N9', "C1'", "C2'"], 'rG': ['N9', "C1'", "C2'"] }
    else:
        av = defaultdict(lambda: on)

    dv = av[res_dir.resname.strip()]
    rv = av[res_ref.resname.strip()]

    dir_points = np.array([res_dir[v].coord for v in dv])
    ref_points = np.array([res_ref[v].coord for v in rv])

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
            fbc.Configuration.template_residues = bpdb.PDBParser().get_structure('t', data_file(conf.Configuration.template_residue_fn))
    s1 = fbc.Configuration.template_residues
    tchain = list(s1.get_chains())[0]


    templates = { 'A': tchain[1], 'C': tchain[2], 'G': tchain[3], 'U': tchain[4]}

    for chain_name, chain in chains.items():
        residues = chain.get_list()

        for residue in residues:
            if residue.id[0].strip():
                residue.id = (" ", residue.id[1], residue.id[2])
            #num = ress[i].id[1]
            old_name = residue.resname.strip()
            target_name = cg.seq.with_missing[fgb.RESID(chain = chain_name, resid = residue.id)]
            log.debug("resname is %s, target name is %s for %s", old_name, target_name, residue)
            # Also replace identical, because they may miss atoms
            #if target_name == old_name:
            #    #Don't replace the correct residues
            #    log.debug("Not replacing %s by %s", old_name, target_name )
            #    continue
            log.debug("Replacing %s by %s", old_name, target_name)
            ref_res = templates[target_name]
            new_res = align_residues(residue, ref_res)
            nsca =  ftup.nonsidechain_atoms + ["OP1", "OP2", "H5'", "H3'", "H4'", "H2'", 'H5"']
            log.debug("NSCA: %s", nsca)
            for atom in residue.child_list[:]: # Need to iterate over copy
                log.debug("%s", atom.get_name() )
                if atom.get_name() not in nsca:
                    log.debug("Detaching %s", atom)
                    residue.detach_child( atom.get_name())
                else:
                    log.debug("Not detaching %s", atom)
                    assert atom.get_name() not in ftup.side_chain_atoms[old_name]
                    assert atom.get_name() not in ftup.side_chain_atoms[target_name]
            sca = ftup.side_chain_atoms[target_name]
            for aname in sca:
                log.debug("Now adding %s", aname)
                residue.add(new_res[aname])

            residue.resname = target_name

            if "O3'" not in residue:
                log.info("Insertion O3' atom to %s", residue)
                new_res = align_residues(residue, ref_res, on=["C3'", "C4'", "C5'"])
                residue.add(new_res["O3'"])



def reorder_residues(chains, cg):
    '''
    Reorder the nucleotides in the chain's list so that they match the order
    in the cg representation.

    :param chains: A dict {chain_id:Bio.PDB.Chain}
    :param cg: A coarse grain representation
    '''
    for chain_name, chain in chains.items():
        chain.child_list.sort(key=lambda x: x.id)
