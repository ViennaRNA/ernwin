#!/usr/bin/python

from optparse import OptionParser

import numpy as np
import math
import random as rand
import sys

import fess.builder.models as cbm
import fess.builder.reconstructor as cbr
import forgi.threedee.model.coarse_grain as ftmc
import forgi.threedee.model.stats as cbs
import forgi.threedee.utilities.graph_pdb as cgg
import forgi.threedee.utilities.vector as cuv
import forgi.utilities.debug as cud
import forgi.utilities.stuff as fus

def get_random_stem_stats():
    '''
    Get a random stem from the statistics and return it.

    @return: A StemModel corresponding to a random stem
    '''
    # pick a length for the stem
    stem_stats = cbs.get_stem_stats()

    l = rand.choice(stem_stats.keys())
    return rand.choice(stem_stats[l])

def get_random_angle_stat(min_len=0., max_len=100.):
    '''
    Create a random angle stats. This refers to the orienation
    of one helix with respect to another.

    @param min_len: The minimum separation between the two stems.
    @param max_len: The maximum separation between the two stems.

    @return: A random AngleStat
    '''
    a =  cbs.AngleStat('', 0, 0,
                             rand.uniform(0, math.pi),
                             rand.uniform(0, 2 * math.pi),
                             rand.uniform(0, 2 * math.pi),
                             rand.uniform(min_len, max_len),
                             rand.uniform(0, math.pi),
                             rand.uniform(0, 2 * math.pi))
    return a

def construct_test_graph(s1_stats, s2_stats, ang_stat, link_length):
    '''
    Construct a simple bulge graph containing two stems, one single
    stranded edge linking them, and a single stranded single connected
    node attached to the first stem.

    @param s1_stats: The statistics describing the first stem
    @param s2_stats: The statistics describing the second stem
    @param link_length: The length of the region linking the two stems

    @return: A BulgeGraph of the described structure
    '''
    s1_len = s1_stats.define[1] - s1_stats.define[0] + 1
    s2_len = s2_stats.define[1] - s2_stats.define[0] + 1

    start_loop_len = 3

    dotbracket = '%s%s%s%s%s%s' % ('.' * start_loop_len, 
                                   '(' * s1_len,
                                   ')' * s1_len,
                                   '.' * link_length,
                                   '(' * s2_len,
                                   ')' * s2_len)
    seq = fus.gen_random_sequence(len(dotbracket))
    cg_db = ftmc.CoarseGrainRNA(dotbracket_str=dotbracket,
                               seq=seq)
    cud.pv('cg_db.to_bg_string()')
    cud.pv('cg_db.get_bulge_dimensions("m0")')

    return cg_db

def bulge_virtual_residue_distance(bg, ld):
    '''
    Calculate the distance between the two virtual residues
    on either side of a bulge.

    @param bg: A BulgeGraph data structure.
    @param ld: The name of the bulge region.

    @return: The distance between the two virtual residues flanking
             flanking the bulge region.
    '''
    if len(bg.edges[ld]) == 2:
        dist2 = cgg.junction_virtual_res_distance(bg, ld)
    else:
        dist2 = 0.

    return dist2

def bulge_virtual_atom_distance(bg, ld):
    '''
    Calculate the distance between the O3' atom and the P
    atom of the two strands that need to be closed.

    @param bg: The BulgeGraph data structure
    @param ld: The name of the bulge region

    @return: The distance between the O3' and P' virtual atoms
             of the flanking residues
    '''
    if len(bg.edges[ld]) == 2:
        dist2 = cgg.junction_virtual_atom_distance(bg, ld)
    else:
        dist2 = 0.

    return dist2

def bulge_length(bg, ld):
    '''
    Calculate the physical length of a bulge region. This is equal to
    the distance between the ends of the two stems that flank the region.
    '''
    connecting_stems = list(bg.edges[ld])

    (s1b, s1e) = bg.get_sides(connecting_stems[0], ld)
    (s2b, s2e) = bg.get_sides(connecting_stems[1], ld)

    return cuv.vec_distance(bg.coords[connecting_stems[0]][s1b],
                             bg.coords[connecting_stems[1]][s2b])

def main():
    parser = OptionParser()

    parser.add_option('-n', '--nucleotides', dest='num_nucleotides', 
                      default=3, help="The number of nucleotides to sample", 
                      type='int')
    parser.add_option('-i', '--iterations', dest='iterations',
                      default=1, help="The number of times to repeat the \
                      simulation", type='int')
    parser.add_option('-a', '--all_lengths', dest='all_lengths', default=False,
                      action='store_true', help='Try all nucleotide lengths up \
                      to the value specified by the --nucleotides parameter.')
    parser.add_option('-o', '--output_file', dest='output_file', default=None,
                    help="The file to dump the output to", type='string')

    (options, args) = parser.parse_args()

    if options.all_lengths:
        start_len = 1
    else:
        start_len = options.num_nucleotides

    if options.output_file != None:
        output = open(options.output_file, 'w')
    else:
        output = sys.stdout

    for k in xrange(start_len, options.num_nucleotides+1):
        min_len = 0. + 3. * k
        max_len = 12. + 7. * k

        for d in np.linspace(min_len, max_len, options.iterations):
            s1 = get_random_stem_stats()
            s2 = get_random_stem_stats()

            ang_stat = get_random_angle_stat(min_len= d, 
                                             max_len = d)

            # Construct a simple graph with an edge of length 3
            bg = construct_test_graph(s1, s2, ang_stat, k)

            # Create a spatial model that will be used to create
            # the 3D model
            try:
                sm = cbm.SpatialModel(bg)
            except IndexError as ie:
                # This can be caused by trying to sample a junction region
                # which is too long and we don't have statistics for
                print >>sys.stderr, "Index error in cbm.SpatialModel(bg)"
                print >>sys.stderr, ie
                continue

            # Indiciate which statistics to use for the 3D model construction
            sm.sample_stats()

            sm.stem_defs['s0'] = s1
            sm.stem_defs['s1'] = s2
            sm.angle_defs['m0'][2] = ang_stat
            sm.angle_defs['m0'][4] = ang_stat
            sm.angle_defs['m0'][5] = ang_stat
            sm.angle_defs['m0'][6] = ang_stat

            # Create the model
            sm.traverse_and_build()
            #sm.bg.output('out.bg')

            # Reconstruct the stems
            try:
                chain = cbr.reconstruct_stems(sm)
            except IOError as ie:
                # Sometimes we'll be missing a fragment
                # Just issue a warning and keep on truckin'
                print >>sys.stderr, "Missing fragment..."
                print >>sys.stderr, ie
            except KeyError as ke:
                print >>sys.stderr, "KeyError in reconstructing stems..."
                print >>sys.stderr, ke

            try:
                ((a,b,i1,i2), best_loop_chain, min_dist) = cbr.reconstruct_loop(chain, sm, 'm0', side=0, samples=3, 
                                     consider_contacts=False, consider_starting_pos=False)

                # Calculate the distances in the coarse grain model
                dist1 = bulge_length(bg, 'm0')
                dist2 = bulge_virtual_residue_distance(bg, 'm0')
                dist3 = bulge_virtual_atom_distance(bg, 'm0')

                output.write("%d %f %f %f %f\n" % (k, dist1, dist2, dist3, min_dist))
                output.flush()
                #print k, dist1, dist2, dist3, min_dist
            except Exception as e:
                print >>sys.stderr, e

if __name__ == "__main__":
    main()

