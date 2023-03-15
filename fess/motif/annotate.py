from __future__ import print_function
from __future__ import absolute_import
import os.path as op
import os
import subprocess as sp
import pandas as pa
import warnings
from . import motif_atlas as ma
import collections as clcs
import fess.builder.config as cbc
import forgi.utilities.debug as fud
import forgi.threedee.model.coarse_grain as ftmc
import sys
import forgi.graph.bulge_graph as fgb
import logging
from six.moves import map

log = logging.getLogger(__name__)
all = [ "annotate_structure" ]

JARED_DIR = op.expanduser(cbc.Configuration.jar3d_dir)
JARED_BIN = cbc.Configuration.jar3d_jar
IL_FILE = cbc.Configuration.jar3d_IL #Download from http://rna.bgsu.edu/data/jar3d/models/ #Relative to JARED_DIR
MOTIF_ATLAS_FILE = cbc.Configuration.jar3d_motif #Click Download at http://rna.bgsu.edu/rna3dhub/motifs/release/il/current# #Relative to JARED_DIR
def annotate_structure(cg, temp_dir, exclude_structure=None, jared_file=None, il_file=None, atlas_file=None):
    '''
    Get the motifs present in this structure.

    :param cg: A CoarseGrainRNA
    :param temp_dir: A directory to place the intermediate files
    :param exclude_structure: None or a string containing a pdb id.
    :param jared_file: path to the jared executable
    :param il_file: path to the interior loop motif atlas file.

    :return: A string containing the motifs.
    '''
    temp_dir = op.expanduser(temp_dir)
    # enumerate the interior loops in the structure
    loop_file = op.join(temp_dir, 'loops')
    try:
        os.makedirs(op.dirname(loop_file))
    except OSError:
        pass
    with open(loop_file, 'w') as f:
        loop_str = cg_to_jared_input(cg)
        f.write(loop_str)

    #fud.pv('jared_file')
    if jared_file is None:
        jared_file = op.expanduser(op.join(JARED_DIR,JARED_BIN))

    if il_file is None:
        il_file = op.expanduser(op.join(JARED_DIR,IL_FILE))

    # run the loops through JAR3D
    jared_output = op.join(temp_dir, 'jared_output')
    cmd = ['java', '-jar', jared_file,
              loop_file, il_file,
              op.join(temp_dir, 'IL_loop_results.txt'),
              op.join(temp_dir, 'IL_sequence_results.txt')]

    #fud.pv("cmd")
    #fud.pv('" ".join(cmd)')
    devnull = open('/dev/null', 'w')
    p = sp.Popen(cmd, stdout=devnull)
    out, err = p.communicate()
    return parse_jared_output(op.join(temp_dir, 'IL_sequence_results.txt'), atlas_file,
                       exclude_structure=exclude_structure, cg=cg)



def get_cg_from_pdb(pdb_file, chains, args, temp_dir=None, cg_filename=None):
    '''
    Get a BulgeGraph from a pdb file.

    :param pdb_file: The filename of the pdb file
    :param chains: The chain ids within the file for which to load the BulgeGraph.
                   If more than one chain is given, they must be connected.
    :param cg_filename: If given, write the cg to this file
    '''
    if temp_dir is not None:
        temp_dir = op.join(temp_dir, 'cg_temp')

    log.info("Creating CG RNA for: %s", pdb_file)
    cg, = ftmc.CoarseGrainRNA.from_pdb(pdb_file, load_chains = chains,
                                       remove_pseudoknots = False,
                                       dissolve_length_one_stems=not args.keep_length_one_stems,
                                       annotation_tool=args.pdb_annotation_tool)

    if cg_filename is not None:
        cg.to_file(cg_filename)
    return cg

def cgdirname_from_args(args):
    if args.pdb_annotation_tool:
        annot_tool=args.pdb_annotation_tool
    else:
        import forgi.config
        c = forgi.config.read_config()
        if "PDB_ANNOTATION_TOOL" in c:
            annot_tool = c["PDB_ANNOTATION_TOOL"]
        else:
            log.warning("No preferred PDB-Annotation-tool set. Inconcistencies due to cached data are possible.")
            annot_tool="?" # In this case, inconsistencies are possible.
    return "cgs_{}_{}".format(int(args.keep_length_one_stems), annot_tool)

def get_coarse_grain_files(struct_name, chains, args, temp_dir=None):
    '''
    Load all connected coarse-grain files for a structure.
    Download the corresponding pdb, if needed.

    :param struct_name: The name of the structure (i.e. '1Y26')
    :param chains: A sequence of chain_ids. If more than one chain_id is given,
                   the chains have to be connected by at least one basepair.
    @return: A forgi.graph.bulge_graph structure describing this chain.
    '''
    CG_DIR = op.join(JARED_DIR, cgdirname_from_args(args))
    PDB_DIR = op.join(JARED_DIR, "pdbs")

    if not op.exists(PDB_DIR):
        os.makedirs(PDB_DIR)

    if not op.exists(CG_DIR):
        os.makedirs(CG_DIR)

    cg_filename = op.join(CG_DIR, struct_name+"_"+"-".join(sorted(chains))+".cg")

    # do we already have the cg representation
    if op.exists(cg_filename):
        return ftmc.CoarseGrainRNA.from_bg_file(cg_filename)
    else:
        pdb_filename = op.join(PDB_DIR, struct_name + ".pdb")

        #do we at least have a pdb file
        if op.exists(pdb_filename):
            return get_cg_from_pdb(pdb_filename, chains,
                                   temp_dir=temp_dir, cg_filename=cg_filename, args=args)
        else:
            log.info ("Downloading pdb for: %s", struct_name)
            import six.moves.urllib.request, six.moves.urllib.error, six.moves.urllib.parse
            response = six.moves.urllib.request.urlopen('http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=%s' % (struct_name))
            html = response.read()


            with open(pdb_filename, 'w') as f:
                f.write(html)
                f.flush()

                return get_cg_from_pdb(pdb_filename, chains, temp_dir=temp_dir,
                                      cg_filename=cg_filename, args=args)


def print_stats_for_motifs(motifs, filename, args, temp_dir=None):
    '''
    Convert all of the motif alignments to coarse-grain element names. This
    requires that the coarse grain representation of the pdb file from
    which the motif comes from be loaded and the element name be determined
    from the nucleotides within the motif.

    :param motifs: A dictionary indexed by an element name. The values are the
                   json motif object from the BGSU motif atlas.
    :param filename: The filename where the stats will be written to.
    :param args: Tha argparse Namespace object. Needed, to use the correct PDB annotation tool.
    '''
    new_motifs = clcs.defaultdict(list)
    i=0
    with open(filename, "w") as file_:
        for key in motifs:
            for motif_entry in motifs[key]:
                log.info(motif_entry)
                for a in motif_entry['alignment']:
                    alignment = ma.MotifAlignment(motif_entry['alignment'][a],
                                            motif_entry['chainbreak'])

                    try:
                        cg = get_coarse_grain_files(alignment.struct,
                                                 temp_dir=temp_dir,
                                                 chains = alignment.chains,
                                                 args=args)
                    except fgb.GraphConstructionError as e:
                        log.warning("Skipping JAR3D entry for {}. Could not "
                                    "construct BulgeGraph because: {}".format(alignment, e))
                        continue

                    elements = set()
                    for r in alignment.residues:
                        log.info(r)
                        elements.add(cg.get_elem(r))

                    loop_elements = set()
                    for e in elements:
                        if e[0] != 's':
                            loop_elements.add(e)

                    try:
                        element_id, = loop_elements
                    except (TypeError, ValueError):
                        log.debug("Skipping JAR3D entry for %s. Elements %s in cg do not match JAR3D.",alignment, elements)
                        continue

                    stats = cg.get_stats(element_id)
                    for stat in stats:
                        i+=1
                        # To ensure unique stat-ids, we use 'j' to identify JAR3D followed by an increasing integer.
                        stat.pdb_name = motif_entry["motif_id"]+"_"+stat.pdb_name+":{}_j{}".format(element_id[0], i)
                        print(stat, file = file_)

def cg_to_jared_input(cg):
    '''
    Take a coarse grain RNA and output all of the loop
    regions within it in a format that JAR3D can understand.

    :param cg: A CoarseGrainRNA structure
    :return: A string containing the interior loops for jared
    '''
    bg = cg
    out_str = ''

    #iterate over the interior loops
    loops = False
    for il in bg.iloop_iterator():
        # get a tuple containing the sequence on each strand
        seqs = bg.get_define_seq_str(il, adjacent=True)
        il_id = ">%s_%s" % (bg.name,
                                  "_".join(map(str, bg.defines[il])))
        out_str += il_id + "\n"
        out_str += "*".join(seqs) + "\n"
        loops = True

    if not loops:
        raise ValueError("No interior loops found in structure")

    return out_str

def parse_jared_output(sequence_results, motif_atlas_file=None, exclude_structure=None, cg=None):
    '''
    Parse the output of the JAR3D file and return all of the motifs.

    :param sequence_results: The sequence results file from JAR3D.
    :param motif_atlas_file: The location of the motif atlas.
    '''
    if motif_atlas_file is None:
        motif_atlas_file = op.join(JARED_DIR, MOTIF_ATLAS_FILE)
    motif_atlas_file = op.expanduser(motif_atlas_file)
    #print ("SEQ", sequence_results)
    data = pa.read_csv(sequence_results)
    atlas = ma.MotifAtlas(motif_atlas_file)
    found_motifs = clcs.defaultdict(list)

    for motif in set(data['identifier']): #In older versions of JAR3D, identifier was sequenceId
        subdata = data[data['identifier'] == motif]
        with warnings.catch_warnings():
            #We do not care if subdata is a view or copy from data.
            #We assign to subdata, but never access the corresponding part of data later on!
            warnings.simplefilter("ignore")
            subdata['score'] = subdata['score'].astype(float)
        subdata = subdata.sort_values(by='score', ascending=False)
        for i, row in subdata.iterrows():
            if not row["passedCutoff"]:
                continue
            motif_id = row['motifId'].split('.')[0]
            motif_entry = atlas.motifs[motif_id]
            res_num = int(motif.split('_')[-1])

            if exclude_structure:
                if atlas.struct_in_motif(motif_id, exclude_structure):
                    # this motif comes from the given structure so we'll exclude it
                    # when reporting the results
                    log.warning("Excluding JAR3D hit %s %s %s, because it is from the input structure.", cg.get_node_from_residue_num(res_num), motif_id, motif_entry['common_name'])
                    continue

            if cg:
                #print '--------------------------------'
                element_name = cg.get_node_from_residue_num(res_num)
                #print element_name, motif, motif_id, motif_entry['common_name']

                if motif_entry['alignment']:
                    '''
                    for a in motif_entry['alignment']:
                        # Print out where this motif comes from
                        print ma.MotifAlignment(motif_entry['alignment'][a],
                                                motif_entry['chainbreak'])
                    '''

                    found_motifs[element_name] += [motif_entry]
            else:
                print ('x', motif, motif_id, motif_entry['common_name'], motif_entry['alignment'])

    return found_motifs
