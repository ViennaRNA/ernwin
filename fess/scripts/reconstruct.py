#!python
from __future__ import print_function

from __future__ import absolute_import
import forgi.utilities.commandline_utils as fuc
import forgi.threedee.utilities.pdb as ftup

import fess.builder.models as fbm
import fess.builder.reconstructor as fbr

import fess.builder.stat_container as fbs
import sys
import collections
import logging
import logging_exceptions
import os.path as op
log = logging.getLogger("reconstruct.py")

def get_parser():
    parser = fuc.get_rna_input_parser("Reconstruct a pdb from a forgi *.coord file produced by ernwin.",
                                      "+", "only_cg", enable_logging=True)
    parser.add_argument("--source-pdb-dir", type=str, required=True,
                        help="A directory where all the pdb files, from which fragments "
                             "can be extracted, are located.")
    parser.add_argument("--source-cg-dir", type=str, required=True,
                        help="A directory where the cg-files corresponding to the pdb files in "
                             "source-pdb-dir are located.")
    parser.add_argument("--interactive", action="store_true")
    parser.add_argument("--server", action="store_true")
    parser.add_argument("--reassign-broken",action="store_true")
    fbs.update_parser(parser)
    return parser

def reconstruct(cg, fn, args, rec):
    log.info("Processing %s", fn)
    stat_source = fbs.from_args(args, cg)
    sm = fbm.SpatialModel(cg)
    sm.load_sampled_elems(stat_source)
    if args.reassign_broken:
        sm.bg.traverse_graph()
        for ml in sm.bg.mloop_iterator():
            if ml not in sm.bg.mst:
                try:
                    del sm.elem_defs[ml]
                except KeyError:
                    pass
        fbm._perml_energy_to_sm(sm, "MAX10000[1FJC1]", stat_source)
        for ml in sm.bg.mloop_iterator():
            if ml not in sm.bg.mst:
                e = sm.junction_constraint_energy[ml].eval_energy(sm.bg)
                log.info("Deviation for %s is %s", ml, e)
                energies = list(sm.junction_constraint_energy[ml].iterate_energies())
                if len(energies)==1:
                    log.debug("Getting used_stat for single energy %s",energies)
                    used_stat = energies[0].used_stat
                else:
                    log.debug("Multiple energies. Iterating")
                    for ef in energies:
                        if ef.element==ml:
                            used_stat = ef.used_stat
                log.info("Used stat for %s is %s", ml, used_stat)
                sm.elem_defs[ml] = used_stat
        sm.save_sampled_elems()
        log.info("Broken stats reassigned. Storing cg.")
        sm.bg.to_file(fn+".reassigned.cg")
    sm.new_traverse_and_build()
    chains = rec.reconstruct(sm)
    print("Writing", fn+".reconstr.pdb")
    chains = ftup.rename_chains_for_pdb(chains)
    ftup.output_multiple_chains(list(chains.values()), fn+".reconstr.pdb")




def main(args):
    rec = fbr.Reconstructor(args.source_pdb_dir, args.source_cg_dir, args.server)
    with fuc.hide_traceback(): # Applies only to WrongFileFormat
        cgs, fns = fuc.cgs_from_args(args, rna_type="only_cg", enable_logging=True, return_filenames=True)
    # Preprocessing
    most_common_pdbs = collections.Counter()
    for cg in cgs:
        sm = fbm.SpatialModel(cg)
        sm.load_sampled_elems(None)
        curr_fns = set()
        for stat in sm.elem_defs.values():
            stat_name = stat.pdb_name
            pdb_basename = stat_name.split(":")[0]
            pdb_filename = op.expanduser(op.join(rec.pdb_library_path, "_".join(pdb_basename.split("_")[:-1])+".cif"))
            try:
                with open(pdb_filename): pass
            except IOError:
                pdb_filename = pdb_filename.rstrip(".cif")+".pdb"
            curr_fns.add(pdb_filename)
        for fn in curr_fns:
            most_common_pdbs[fn]+=1
    for fn, count in most_common_pdbs.most_common(250):
        if count==1:
            break
        print("Preloading {}, used {} times".format(fn, count))
        rec.get_pdb(fn, True)
    print("Preloading of most common PDBs done")
    logging.getLogger("forgi").setLevel(logging.ERROR)
    logging_exceptions.config_from_args(args)
    for i, cg in enumerate(cgs):
        try:
            fn = fns[i]
            reconstruct(cg, fn, args, rec)
        except Exception as e:
            logging_exceptions.log_exception(e)
            log.exception("During reconstruction of cg %s, an error occurred: %s", fn, e)
    i=0
    #print("Waiting for stdin", file=sys.stderr)
    #while True:
    #    try:
    #        line = sys.stdin.readline()
    #        fn=line.strip()
    #        print("Received", fn, file=sys.stderr)
    #        cg = fuc.load_rna(fn, "only_cg", False )
    #        reconstruct(cg, fn, args, rec)
    #    except Exception as e:
    #        logging_exceptions.log_exception(e)
    #        log.exception("During reconstruction of additional cg %s, an error occurred: %s", fn, e)



parser = get_parser()
if __name__=="__main__":
    args = parser.parse_args()
    main(args)
