from __future__ import print_function
from __future__ import absolute_import
import argparse
import sys
from collections import defaultdict
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.utilities.commandline_utils as fuc
import logging

def get_parser():
    parser=fuc.get_rna_input_parser("Create stats", nargs='+', rna_type="only_cg")
    return parser
parser=get_parser()

if __name__ == "__main__":
    next_id = defaultdict(int)
    args = parser.parse_args()
    cgs = fuc.cgs_from_args(args, "only_cg")
    for cg in cgs:
        if sys.stderr.isatty():
            print(cg.name, file=sys.stderr)
        for elem in cg.defines.keys():
            if elem in cg.incomplete_elements:
                print("Skipping element", elem, file = sys.stderr)
                continue
            base_name = "{}:{}_".format(cg.name, elem[0])
            for stat in cg.get_stats(elem):
                if len(stat.vbase)  == len(stat.vsugar) == len(stat.vbackbone) == len(list(cg.define_residue_num_iterator(elem))):
                    idnr = next_id[base_name]
                    next_id[base_name]+=1
                    name = base_name+str(idnr)
                    stat.pdb_name = name
                    if elem.startswith("m"):
                        try:
                            dist = ftug.junction_virtual_atom_distance(cg, elem)
                            stat_dist = stat.get_virtual_atom_distance()
                        except:
                            print(stat)
                            #raise
                        else:
                            print(stat, "# distance: {}. stat_dist {}".format(dist, stat_dist))
                    else:
                        print(stat)
                else:
                    print("#IGNORE:", stat, "3 lengths", len(stat.vbase),
                          len(stat.vsugar), len(stat.vbackbone),
                          len(list(cg.define_residue_num_iterator(elem))),  file = sys.stderr)
