from __future__ import print_function
import argparse
import sys
from collections import defaultdict
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.pdb as ftup
import forgi.graph.residue as fgr
import forgi.utilities.commandline_utils as fuc
import os.path as op
import logging

def generateParser():
    parser=fuc.get_rna_input_parser("Extract fragments",
                                     nargs="+", rna_type="pdb", enable_logging=True)
    return parser
parser=generateParser()

if __name__ == "__main__":
    next_id = defaultdict(int)
    args = parser.parse_args()
    cgs, filenames = fuc.cgs_from_args(args, "+", "pdb", enable_logging=True,
                                       return_filenames=True, skip_errors=True)

    for cg, fn in zip(cgs, filenames):
        if sys.stderr.isatty():
            print(cg.name, file=sys.stderr)
        for elem in cg.defines:
            if elem in cg.incomplete_elements:
                continue
            base_name = "{}:{}_".format(cg.name, elem[0])
            idnr = next_id[base_name]
            next_id[base_name]+=1
            name = base_name+str(idnr)
            fragment_chains = ftup.extract_subchains_from_seq_ids(cg.chains,
                                    cg.define_residue_num_iterator(elem, seq_ids=True, adjacent=elem[0]!="s"))
            try:
                ftup.output_multiple_chains(list(fragment_chains.values()), op.join("fragments", name+".pdb"))
            except AttributeError:
                continue
            direction = "a"
            remark=""
            for resid in cg.define_residue_num_iterator(elem, seq_ids=True, adjacent=elem[0]!="s"):
                if resid in cg.seq._modifications:
                    remark+="# {}={}".format(fgr.resid_to_str(resid), cg.seq._modifications[resid])
            for stat in cg.get_stats(elem):
                stat_name = name+direction
                stat.pdb_name = stat_name
                print(stat, remark)
                direction="b"
