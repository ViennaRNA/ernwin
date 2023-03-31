#!python
from __future__ import absolute_import, division, print_function, unicode_literals
from builtins import (ascii, bytes, chr, dict, filter, hex, input, #pip install future
                      int, map, next, oct, open, pow, range, round,
                      str, super, zip) #future package

import argparse
import sys
import random
import os
import copy
import multiprocessing
import traceback
import itertools
import numpy as np

import forgi.utilities.commandline_utils as fuc
import forgi.threedee.utilities.vector as ftuv
import fess.builder.stat_container as fbstat
import fess.directory_utils
import fess.builder.energy as fbe
from fess.builder import config
import fess.builder.monitor as fbm
import fess.builder.move as fbmov
import fess.builder._other_movers
import fess.builder.replicaExchange as fbr
import fess.builder.builder as fbb
import fess.builder.models as fbmodel
import fess.builder.sampling as fbs
from fess.utils import get_version_string

import logging
from six.moves import range

log=logging.getLogger("ernwin.__main__")

parser = fuc.get_rna_input_parser("ERNWIN: Coarse-grained sampling of RNA 3D structures.",
                                  nargs=1, rna_type="cg",
                                  parser_kwargs={"formatter_class":argparse.RawTextHelpFormatter})
parser.add_argument('--seed', action='store',
                    help="Seed for the random number generator. It is taken module 2**32 to work with numpy.",
                    type=int)
parser.add_argument('-i', '--iterations', action='store', default=10000, help='Number of structures to generate', type=int)
parser.add_argument('--replica-exchange', type=int, help="Experimental")
parser.add_argument('--replica-exchange-every-n', type=int, default=1, help="After how many steps try a replica exchange?")
parser.add_argument('--parallel', action="store_true",
                    help="Spawn parallel processes.\n"
                         "Only used if either --replica-exchange or --num_builds\n"
                         "is given. Perform sampling in parallel\n"
                         "with one process per replica/ build.")
parser.add_argument('--num-builds', action='store', type=int,
                    default=1,
                    help="Number of RNA structures built independently.\n"
                         "Each of it will be sampled for \n"
                         "--iterations steps independently")

parser.add_argument('--externally-interacting', type=str,
                            help="A comma-separated list of element names or nt positions, "
                                 "which have interactions with proteins/... and "
                                 " will be excluded from A-Minor and "
                                 "loop-loop interaction energies.")

parser.add_argument('--reference', action='store', type=str,
                    help="A *.cg/ *.coord or *.pdb file.\n"
                         "Calculate the RMSD and MCC relative to the structure\n"
                         "in this file, not to the structure used as starting\n"
                         "point for sampling.")


# Each of the following modules of ernwin adds its own options to the parser.
for module in [fbstat, fess.directory_utils, fbe, fbmov, fbm, fbb, fbmodel]:
    module.update_parser(parser)


def main(args):
    args = parser.parse_args(args)

    cg, = fuc.cgs_from_args(args, rna_type="cg",
                            enable_logging=True) # Set log-level as a sideeffect

    if args.reference:
        args.rna = [args.reference]
        reference_cg, = fuc.cgs_from_args(args, rna_type="cg", enable_logging=False)
    else:
        reference_cg = None
    if args.externally_interacting:
        elems_or_nts=args.externally_interacting.split(",")
        nts = []
        for pos in elems_or_nts:
            try:
                nt = int(pos)
            except Exception:
                d = cg.defines[pos] #Might raise a KeyError, but that is fine
                if d:
                    nts.append(d[0])
                else:
                    log.warning("Not setting %s to externally interacting, because it contains 0 nucleotides.", elem)
            else:
                nts.append(nt)
        cg.interacting_residues.extend([ cg.seq.to_resid(nt) for nt in nts])
        log.info("The RNA now has the following elements not perticipating in interaction energies (if presenty): %s", cg.interacting_elements)

    if len(list(cg.stem_iterator()))<2:
        raise ValueError("No sampling can be done for structures with fewer than 2 stems")

    with fess.directory_utils.make_outdir(args, cg) as main_dir:
        cg_stri = cg.to_cg_string()
        with open(os.path.join(main_dir, 'input.cg'), "w") as f:
            print(cg_stri, file=f)
        try:
            run(args, cg, main_dir, reference_cg)
        except BaseException as e:
            with open(os.path.join(main_dir, 'exception.log'), "w") as f:
                print("Running on python {}, the following error occurred:".format(sys.version), file=f)
                print("{}: {}".format(type(e).__name__, str(e)), file=f)
                print(str(traceback.format_exc()), file=f)
            raise

def run(args, cg, main_dir, reference_cg):
    setup_rng(args)
    stat_source = fbstat.from_args(args, cg) #Uses sampling_output_dir

    # If we perform normal sampling with parallel=True,
    # we start sampling while we are building.
    # Otherwise (replica-exchange, single-process sampling),
    # we sample after all structures were built.
    samplers = []
    processes = []
    if reference_cg is None:
        reference_cg = cg
    for i, sm in enumerate(build_spatial_models(args, cg, stat_source, main_dir)):
        cg_stri = sm.bg.to_cg_string()
        with open(os.path.join(main_dir,
                               'build{:06d}.coord'.format(i+1)), "w") as f:
            print(cg_stri, file=f)
        # The PDD energy may want to use the original cg
        sampling_energy = fbe.from_args(args, sm.bg, stat_source, i, reference_cg)
        print("Original energy: ", sampling_energy.eval_energy(sm.bg))
        if args.iterations>0:
            sampler = setup_sampler(args, sm, stat_source, sampling_energy, i, reference_cg)
            samplers.append(sampler)
            if not args.replica_exchange and args.parallel:
                p = multiprocessing.Process(target=sample_one_trajectory,
                                            args = (sampler, args.iterations))
                processes.append(p)
                p.start()
    if not args.parallel and not args.replica_exchange:
        for sampler in samplers:
            sample_one_trajectory(sampler, args.iterations)
    elif args.replica_exchange:
        if args.parallel:
            fbr.start_parallel_replica_exchange(samplers, args, stat_source)
        else:
            re = fbr.ReplicaExchange(samplers)
            re.run(args.iterations)

    # Join all processes, if any were started.
    for p in processes:
        p.join()

#if args.reconstruct:
#    raise NotImplementedError("TODO")

def sample_one_trajectory(sampler, iterations):
    with sampler.stats_collector.open_outfile():
        for i in range(iterations):
            sampler.step()
        sampler.stats_collector.collector.to_file()

def build_spatial_models(args, cg, stat_source, main_dir):
    """
    An iterator over spatial models, which contain deepcopies of cg.
    """
    if args.replica_exchange and args.num_builds>1:
        raise ValueError("--replica-exchange and --num-builds are mutually exclusive.")
    if args.replica_exchange:
        build_count = args.replica_exchange
    else:
        build_count = args.num_builds
    build_function = fbb.from_args(args, stat_source, main_dir)
    sm = None
    for i in range(build_count):
        curr_cg = copy.deepcopy(cg)
        sm = fbmodel.from_args(args, curr_cg, stat_source, i)
        build_function(sm)
        yield sm



def setup_sampler(args, sm, stat_source, sampling_energy, replica_nr=None, original_cg=None):
    """
    Part of the setup, that has to be repeated for every
    replica in replica-exchange Monte Carlo.

    :param replica_nr: None, if no replica-exchange is desired.
                       An integer (starting with 0) for the number
                       of the replica.
    """
    # The monitor uses the original structure as reference, IF it has 3D coordinates
    if original_cg.coords.is_filled and original_cg.twists.is_filled:
        show_min_rmsd=True
    else:
        show_min_rmsd=False
    mover = fbmov.from_args(args, stat_source, sm, replica_nr)
    if args.replica_exchange:
        out_dir = os.path.join(config.Configuration.sampling_output_dir,
                          "temperature_{:02d}".format(replica_nr+1))
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
    else:
        out_dir = os.path.join(config.Configuration.sampling_output_dir,
                          "simulation_{:02d}".format(replica_nr+1))
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

    monitor = fbm.from_args(args, original_cg, sampling_energy, stat_source, out_dir, show_min_rmsd)
    sampler = fbs.MCMCSampler(sm, sampling_energy, mover, monitor)
    return sampler

def setup_rng(args):
    if args.seed:
        seed_num=args.seed % 2**32
    else:
        # We need a seed that we can use and print to the user, so the user can reproduce this run in the future.
        seed_num = random.randint(0,4294967295) #4294967295 is maximal value for numpy
    random.seed(seed_num)
    np.random.seed(seed_num)
    with open(os.path.join(config.Configuration.sampling_output_dir, "version.txt"), "w") as versionfile:
        print ("Random Seed: {}".format(seed_num), file=versionfile)
        label = get_version_string()
        print ("Version: {}".format(label), file=versionfile)
        print ("Command:\n  {}".format(" ".join(sys.argv)), file=versionfile)


if __name__ == "__main__":
    main(args=sys.argv[1:])
