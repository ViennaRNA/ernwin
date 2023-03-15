from __future__ import absolute_import
import contextlib
import os.path
import logging

from forgi.utilities.commandline_utils import open_for_out


from .builder import config

log = logging.getLogger(__name__)


def update_parser(parser):
    """
    Add arguments related to output files and directories to the parser class.

    :param parser: A argparse.ArgumentParser instance.
    """
    output_options = parser.add_argument_group(
                                    "Controlling output",
                                    description="These options control the output or ERNWIN")
    #Controll output files
    output_options.add_argument('--output-base-dir', action='store', type=str, default="",
                                help="Base dir for the output. \n"
                                      "In this directory, a subfolder with the name\n"
                                      "from the fasta-file will be generated")
    output_options.add_argument('-o', '--output-dir-suffix', action='store', type=str,
                                default="",
                                help="Suffix attached to the name from the fasta-file, \n"
                                     "used as name for the subfolder with all structures.")



@contextlib.contextmanager
def make_outdir(args, cg):
    """
    A context manager that creates the required output directories
    from args and yields the output file (inside those) opened for writing.

    :param args: A agrparse.Namespace instance, as returned by calling parse_args() on an
                  ArgumentParser instance
    """
    if args.output_base_dir:
        if not os.path.exists(args.output_base_dir):
            os.makedirs(args.output_base_dir)
            log.info("Directory %s created.",args.output_base_dir)
        else:
            log.warning("WARNING: Using existing directory %s. "
                        "Potentially overwriting its content.", args.output_base_dir)
    subdir=cg.name+args.output_dir_suffix
    old_config_dir = config.Configuration.sampling_output_dir
    try:
        config.Configuration.sampling_output_dir = os.path.join(args.output_base_dir, subdir)
        if not os.path.exists(config.Configuration.sampling_output_dir):
            os.makedirs(config.Configuration.sampling_output_dir)
            log.info("Directory %s created. This folder will be used for all output "
                     "files.", config.Configuration.sampling_output_dir)
        else:
            log.warning("Using existing directory %s. "
                        "Potentially overwriting its content.",
                        config.Configuration.sampling_output_dir)
        yield config.Configuration.sampling_output_dir
    finally:
        config.Configuration.sampling_output_dir = old_config_dir
