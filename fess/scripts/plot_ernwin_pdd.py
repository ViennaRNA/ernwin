#!python
from __future__ import division
from __future__ import absolute_import
from __future__ import print_function
import sys
import matplotlib.pyplot as plt
import numpy as np

import fess.builder.energy as fbe
import forgi.utilities.commandline_utils as fuc
import forgi

import argparse
from six.moves import range


def get_pdd(energy, rna):
    m = energy.get_pdd(rna, energy._level, energy._stepwidth, energy.only_seqids)
    print(m)
    m=m[1]
    m = energy.pad(m)
    return m/np.sum(m)
            
if __name__ == "__main__":
    parser = fuc.get_rna_input_parser("Plot pair distance distributions.",
                                      nargs="+", rna_type="cg")
    fbe.update_parser(parser)     
    
    args = parser.parse_args()                                 
    rnas = fuc.cgs_from_args(args, rna_type="cg",
                              enable_logging=True) # Set log-level as a sideeffect
    combined_energy = fbe.from_args(args, rnas[0], None, None, None)
    
    pdd_energies = [energy for energy in combined_energy.iterate_energies()
                    if isinstance(energy, (fbe.PDDEnergy, fbe.Ensemble_PDD_Energy))]

    for energy in pdd_energies:
        plt.title(energy.shortname)
        plt.plot([energy._stepwidth*i for i in range(len(energy.target_values))], 
                  energy.target_values, label="Experiment", color="red")
        if isinstance(energy, fbe.Ensemble_PDD_Energy):
            error = energy.check_target_pdd(energy.target_pdd["distance"], 
                                            energy.error, normalize=False)
            plt.plot([energy._stepwidth*i for i in range(len(energy.target_values))], 
                      energy.target_values + error/energy.normalization_factor*energy.adjustment, 
                      color="orange", label="Target Standard deviation (adjustment)")
            plt.plot([energy._stepwidth*i for i in range(len(energy.target_values))],
                      np.maximum(0, energy.target_values - error/energy.normalization_factor*energy.adjustment), 
                      color="orange")
                      
        for rna_num, rna in enumerate(rnas):
            plt.plot([energy._stepwidth*i for i in range(len(energy.target_values))], 
                     get_pdd(energy, rna), "o", label="predicted structure  {}".format(rna_num))
        plt.xlabel("r [A]")
        plt.legend()
        plt.show()
