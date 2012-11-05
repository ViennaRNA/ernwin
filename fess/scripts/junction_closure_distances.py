#!/usr/bin/python

import sys
import random as rand
import math
from optparse import OptionParser

import corgy.builder.models as cbm
import corgy.builder.reconstructor as cbr
import corgy.builder.stats as cbs
import corgy.graph.bulge_graph as cgb
import corgy.graph.graph_pdb as cgg
import corgy.utilities.debug as cud
import corgy.utilities.vector as cuv

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
    return cbs.AngleStat('', 0, 0,
                             rand.uniform(0, math.pi),
                             rand.uniform(0, 2 * math.pi),
                             rand.uniform(0, 2 * math.pi),
                             rand.uniform(min_len, max_len),
                             rand.uniform(0, math.pi),
                             rand.uniform(0, 2 * math.pi))

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
    s1_len = s1_stats.define[1] - s1_stats.define[0]
    s2_len = s2_stats.define[1] - s2_stats.define[0]

    # the new defines for the stem
    new_d1 = range(4)
    new_d2 = range(4)

    #Start with a little loop so that the traversal can get started
    start_loop_len = 3
    loop_d = [0, start_loop_len]

    # shift the define of the first stem so that it is at the beginning
    new_d1[0] = 3
    new_d1[1] = 3 + s1_len

    # no residues between the two strands of stem
    new_d1[2] = new_d1[1] + 1
    new_d1[3] = new_d1[2] + s1_len

    # add the junction between the two stems
    new_d2[0] = new_d1[3] + link_length + 1
    new_d2[1] = new_d2[0] + s2_len

    new_d2[2] = new_d2[1] + 1
    new_d2[3] = new_d2[2] + s2_len


    # create the bulge graph programmatically
    bg = cgb.BulgeGraph()
    bg.seq = "".join([rand.choice(['A','C','G','U']) 
                      for i in xrange(new_d2[3])])

    bg.defines['s1'] = new_d1
    bg.defines['s2'] = new_d2
    bg.defines['b1'] = [new_d1[3], new_d2[0]]
    bg.defines['b2'] = loop_d

    bg.edges['s1'] = set(['b1', 'b2'])
    bg.edges['s2'] = set(['b1'])
    bg.edges['b1'] = set(['s1', 's2'])
    bg.edges['b2'] = set(['s1'])

    bg.weights['s1'] = 2
    bg.weights['s2'] = 2
    bg.weights['b1'] = 1
    bg.weights['b2'] = 1

    return bg

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


    (options, args) = parser.parse_args()

    if options.all_lengths:
        start_len = 1
    else:
        start_len = options.num_nucleotides

    for k in xrange(start_len, options.num_nucleotides+1):
        for i in xrange(options.iterations):
            s1 = get_random_stem_stats()
            s2 = get_random_stem_stats()

            ang_stat = get_random_angle_stat(min_len=0., 
                                             max_len = 20. + 15. * k)

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
            sm.stem_defs['s1'] = s1
            sm.stem_defs['s2'] = s2
            sm.angle_defs['b1'] = ang_stat

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
                continue
            
            c, min_dist = cbr.reconstruct_loop(chain, sm, 'b1', side=0, samples=40, 
                                 consider_contacts=False)

            # Calculate the distances in the coarse grain model
            dist1 = bulge_length(bg, 'b1')
            dist2 = bulge_virtual_residue_distance(bg, 'b1')

            print k, dist1, dist2, min_dist

if __name__ == "__main__":
    main()

