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
from contextlib import nested

import numpy as np

import forgi.utilities.commandline_utils as fuc

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

parser = fuc.get_rna_input_parser("MULTIERNWIN: Coarse-grained sampling of RNA 3D structures from multiple secondary structures.",
                                  nargs='+', rna_type="cg",
                                  parser_kwargs={"formatter_class":argparse.RawTextHelpFormatter, "conflict_handler": "resolve"})
parser.add_argument('--seed', action='store',
                    help="Seed for the random number generator. It is taken module 2**32 to work with numpy.",
                    type=int)
parser.add_argument('-i', '--iterations', action='store', default=10000, help='Number of structures to generate', type=int)
parser.add_argument('--shared-energy', action='store', type=str, help='Energy function that is shared between all samplers')
# Each of the following modules of ernwin adds its own options to the parser.
for module in [fbstat, fess.directory_utils, fbe, fbmov, fbm, fbb, fbmodel]:
    module.update_parser(parser)
#Remove some arguments
parser.add_argument('--jar3d', help=argparse.SUPPRESS)


def main():
    args = parser.parse_args()

    cgs = fuc.cgs_from_args(args, rna_type="cg",
                            enable_logging=True) # Set log-level as a sideeffect
    for cg in cgs:
        if len(list(cg.stem_iterator()))<2:
            raise ValueError("No sampling can be done for structures with fewer than 2 stems")

    with fess.directory_utils.make_outdir(args, cgs[0]) as main_dir:
        for i,cg in enumerate(cgs):
            cg_stri = cg.to_cg_string()
            with open(os.path.join(main_dir, 'input{}.cg'.format(i)), "w") as f:
                print(cg_stri, file=f)
        try:
            run(args, cgs, main_dir)
        except BaseException as e:
            with open(os.path.join(main_dir, 'exception.log'), "w") as f:
                print("Running on python {}, the following error occurred:".format(sys.version), file=f)
                print("{}: {}".format(type(e).__name__, str(e)), file=f)
                print(str(traceback.format_exc()), file=f)
            raise

def run(args, cgs, main_dir):
    setup_rng(args)
    stat_source = fbstat.from_args(args, None) #Uses sampling_output_dir

    samplers = []
    shared_energy = None
    for i, cg in enumerate(cgs):
        sm = build_spatial_model(args, cg, stat_source, main_dir)
        cg_stri = sm.bg.to_cg_string()
        with open(os.path.join(main_dir,
                               'build{:06d}.coord'.format(i+1)), "w") as f:
            print(cg_stri, file=f)
        sampler, shared_energy = setup_sampler(args, sm, stat_source, i, shared_energy)
        samplers.append(sampler)

    with nested(*[sampler.stats_collector.open_outfile() for sampler in samplers]):
        for i in range(args.iterations):
            for sampler in samplers:
                sampler.step()


def build_spatial_model(args, cg, stat_source, main_dir):
    """
    An iterator over spatial models, which contain deepcopies of cg.
    """
    build_function = fbb.from_args(args, stat_source, main_dir)
    sm = None
    curr_cg = copy.deepcopy(cg)
    sm = fbmodel.from_args(args, curr_cg, stat_source)
    build_function(sm)
    return sm


def setup_sampler(args, sm, stat_source, replica_nr=None, shared_energies=None):
    """
    Part of the setup, that has to be repeated for every
    replica in replica-exchange Monte Carlo.

    :param replica_nr: None, if no replica-exchange is desired.
                       An integer (starting with 0) for the number
                       of the replica.
    """
    cg=sm.bg
    show_min_rmsd=False
    if shared_energies is None:
        shared_energies=fbe.EnergyFunction.from_string(args.shared_energy,
                                              cg=cg,
                                              stat_source=stat_source,
                                              pdd_target=args.pdd_file,
                                              pdd_stepsize=args.pdd_stepsize)
    sampling_energy = fbe.from_args(args, cg, stat_source, replica_nr)
    sampling_energy = fbe.CombinedEnergy([sampling_energy]+shared_energies)
    mover = fbmov.from_args(args, stat_source, sm, replica_nr)
    out_dir = os.path.join(config.Configuration.sampling_output_dir,
                          "linkedsimulation_{:02d}".format(replica_nr+1))
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    monitor = fbm.from_args(args, cg, sampling_energy, stat_source, out_dir, show_min_rmsd)
    sampler = fbs.MCMCSampler(sm, sampling_energy, mover, monitor, rerun_prev_energy=True)
    return sampler, shared_energies

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
    main()
