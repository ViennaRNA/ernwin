import fess.builder.energy as fbe
import fess.utils as fau
import matplotlib.pyplot as plt
import numpy as np
import math
import logging
import pandas as pd
import argparse
from collections import defaultdict
import forgi.threedee.model.coarse_grain as ftmc
import forgi.utilities.commandline_utils as fuc

fr3d_query_string="""
      pdb = 1ffk
      cutoff = 0.5;
      First nt constrained to A
      2nd-3rd nt constrained to cWW
      NTs (0): 104 (chain 9) 957 1009
      NTs (I): 521 1364 637
      NTs (II): 520 1363 638
      NTs (III): 519 1362 639
      """

def get_parser():
    """
    Here all commandline & help-messages arguments are defined.

    :returns: an instance of argparse.ArgumentParser
    """
    parser = fuc.get_rna_input_parser("Generate target distributions for ernwin", nargs="+")
    ### Options ###
    # Options to modify general behavior
    parser.add_argument('--plot-only', action="store_true", help='Do not generate any files, only plot results from existing files')
    parser.add_argument('--use-subgraphs', type=int, help='Use subgraphs for the target distribution. Default: Use enough subgraphs to generate a total of 10000 graphs.')
    # Files that will be generated
    parser.add_argument('--rog-target-file', help='Filename where the ROG target distribution will be stored', type=str)
    parser.add_argument('--ame-target-file', help='Filename where the AMinor target distribution will be stored', type=str)
    parser.add_argument('--sld-target-file', help='Filename where the Shortest loop distance target distribution will be stored', type=str)
    parser.add_argument('--ame-orientation-outfile', default="stats/test_aminor_orientation_1S72.txt", help='Filename where the AMinor orientations will be stored', type=str)
    parser.add_argument('--chain-id-mapping-dir', help='Directory name where chain-id-mappings for pdb bundles are located.', type=str)
    # Files that are required
    parser.add_argument('--fr3d-result-file', default="fess/stats/AMinor_FR3D_hits.txt", help='Filename with the output of FR3D', type=str)
    #Additional info
    parser.add_argument('--fr3d-query-string', default=fr3d_query_string, help='FR3D query information that will be added to the file generated.', type=str)
    #parser.add_argument('--ame-pdb-id-file', help='For the AMinor energy, only consider pdb ids from this file', type=str)
    parser.add_argument('--precalculated-ame-orient', default=False, action="store_true")
    return parser


parser = get_parser()
if __name__=="__main__":
    args = parser.parse_args()

    #Logging
    logging.basicConfig(format="%(levelname)s:%(name)s:%(funcName)s[%(lineno)d]: %(message)s")
    logging.captureWarnings(True)
    cgs = fuc.cgs_from_args(args, "+")

    log = logging.getLogger(__name__)
    log.info(fau.get_version_string())

    if args.use_subgraphs is None:
        use_subgraphs = max(1, math.ceil(10000/len(cg_files)))
    else:
        use_subgraphs = args.use_subgraphs

    if not args.plot_only:
      #Generating the files

      ## ROG
      if args.rog_target_file:
          fbe.RadiusOfGyrationEnergy.generate_target_distribution(cg_files,
                                              args.rog_target_file,
                                              use_subgraphs = use_subgraphs)

      ## AME
      if args.ame_target_file:
          if args.precalculated_ame_orient:
                all_cgs = defaultdict(list)
                for cg in cgs:
                    all_cgs[cg.name[:4]].append(cg)
                fbe.AMinorEnergy._generate_target_dist_given_orientation(all_cgs, 'fess/'+args.ame_target_file, args.ame_orientation_outfile, use_subgraphs = use_subgraphs)
          else:
              fbe.AMinorEnergy.generate_target_distribution(cgs,
                             args.fr3d_result_file,
                             args.chain_id_mapping_dir,                             
                             out_filename = args.ame_target_file,
                             orientation_outfile = args.ame_orientation_outfile,
                             fr3d_query = args.fr3d_query_string, use_subgraphs = use_subgraphs)

      ## SLD
      if args.sld_target_file:
          fbe.ShortestLoopDistancePerLoop.generate_target_distribution(cgs,
                                                            args.sld_target_file,
                                                            use_subgraphs = use_subgraphs)

    # Do some plotting
    for energy, new_fn in [(fbe.RadiusOfGyrationEnergy, args.rog_target_file),
                           (fbe.ShortestLoopDistancePerLoop, args.sld_target_file)]:
        fig, ax = plt.subplots()
        old = np.genfromtxt(fbe.load_local_data(energy.real_stats_fn), delimiter=' ')
        new = np.genfromtxt(fbe.load_local_data(new_fn), delimiter = ' ')
        try:
          ax.scatter(old[:,0], old[:,1], label = "old")
        except:
          ax.scatter(list(range(len(old))), old, label = "old")
        ax.scatter(new[:,0], new[:,1], label = "new", color = "red")
        ax.legend()
        ax.set_title(energy.__name__)
        plt.show(block=False)

    fig, ax = plt.subplots()
    old = np.genfromtxt(fbe.load_local_data(fbe.AMinorEnergy.real_stats_fn), delimiter=' ')
    new = pd.read_csv(fbe.load_local_data(args.ame_target_file), delimiter=' ', comment = "#")
    print(new)
    print(new.columns)
    print(new[u"rna_length"].dtype)
    new = new.as_matrix([u"rna_length", u"total_prob"])
    print(new)
    ax.scatter(old[:,0], old[:,2], label = "old")
    ax.scatter(new[:,0], new[:,1], label = "new", color="red")
    ax.legend()
    ax.set_title("AME")
    plt.show()
