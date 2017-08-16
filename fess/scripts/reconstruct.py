import forgi.utilities.commandline_utils as fuc
import forgi.threedee.utilities.pdb as ftup

import fess.builder.models as fbm
import fess.builder.reconstructor as fbr

def get_parser():
    parser = fuc.get_rna_input_parser("Reconstruct a pdb from a forgi *.coord file produced by ernwin.",
                                      1, "only_cg", enable_logging=True)
    parser.add_argument("--source-pdb-dir", type=str, required=True,
                        help="A directory where all the pdb files, from which fragments "
                             "can be extracted, are located.")
    parser.add_argument("--source-cg-dir", type=str, required=True,
                        help="A directory where the cg-files corresponding to the pdb files in "
                             "source-pdb-dir are located.")
    parser.add_argument("--outfilename", type=str,
                        help="Target filename for the pdb to be written.", default="out.pdb")
    return parser

def main(args):
    with fuc.hide_traceback(): # Applies only to WrongFileFormat
        cg, = fuc.cgs_from_args(args, 1, rna_type="only_cg", enable_logging=True)
        sm = fbm.SpatialModel(cg)
        sm.load_sampled_elems()
        sm.new_traverse_and_build()
        rec = fbr.Reconstructor(args.source_pdb_dir, args.source_cg_dir)
        chains = rec.reconstruct(sm)
        ftup.output_multiple_chains(chains.values(), args.outfilename)

parser = get_parser()
if __name__=="__main__":
    args = parser.parse_args()
    main(args)
