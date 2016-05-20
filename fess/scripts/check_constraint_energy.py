from __future__ import print_function, absolute_import, division
from builtins import (ascii, bytes, chr, dict, filter, hex, input,
                      map, next, oct, open, pow, range, round,
                      str, super, zip) 


import argparse, random
import forgi.utilities.stuff as fus
import fess.builder.models as fbm
import fess.builder.energy as fbe
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.stats as ftms
import forgi.threedee.utilities.graph_pdb as ftug
from fess import data_file
def generateParser():
    parser=argparse.ArgumentParser( description="For all stats from which we can sample, "
                                                "make sure they fulfill the energy")  
    return parser


def construct_test_graph(s1_len, s2_len, link_length):
    '''
    Construct a simple bulge graph containing two stems, one single
    stranded edge linking them, and a single stranded single connected
    node attached to the first stem.

    @param s1_len: The length of the first stem
    @param s2_stats: The length of the second stem
    @param link_length: The length of the single stranded region linking the two stems

    @return: A BulgeGraph of the described structure
    '''
    start_loop_len = 3

    dotbracket = "".join(['.' * start_loop_len,
                          '(' * s1_len,
                          '.' * 3,
                          ')' * s1_len,
                          '.' * link_length,
                          '(' * s2_len,
                          '.' * 3,
                          ')' * s2_len])
    seq = fus.gen_random_sequence(len(dotbracket))
    cg = ftmc.CoarseGrainRNA(dotbracket_str=dotbracket,
                               seq=seq)
    return cg



parser=generateParser()
if __name__=="__main__":
    args = parser.parse_args()
    energy = fbe.CombinedEnergy([fbe.RoughJunctionClosureEnergy()])
    for link_length in range(0,4):                
        cutoff_distance = (link_length) * 6.22 + 14.0            
        print("\nLink Length {}".format(link_length))
        for l in range(1, 30): #Sampling 30 times with different stem stats
            cg = construct_test_graph(random.randrange(3,15),l,link_length)
            sm=fbm.SpatialModel(cg, ftms.get_conformation_stats(data_file("stats/all_filtered.stats")))
            sm.sample_stats()
            possible_sstats=sm.conf_stats.sample_stats(cg, "s1")
            for sstat in possible_sstats:
                if sstat.pdb_name.startswith("ideal"):
                    break
            maxdist=0
            sm.elem_defs["s1"] = sstat            
            possible_stats=sm.conf_stats.sample_stats(cg, "m0")                
            for i, stat in enumerate(possible_stats):

                #if i%1000 == 0:
                #    print(i, "/", len(possible_stats), end="")
                sm.elem_defs["m0"] = stat
                sm.traverse_and_build()
                e=energy.eval_energy(sm)

                dist = ftug.junction_virtual_atom_distance(cg, "m0")
                if dist>maxdist:
                    maxdist=dist
                if e>0:
                    print ("\nSampled stats do not match energy:")
                    print (stat)
                    print( dist, ">", cutoff_distance)
                    print(sm.elem_defs["s0"])
                    print(sm.elem_defs["s1"])
            print( "Maxdist", maxdist, "<", cutoff_distance)




