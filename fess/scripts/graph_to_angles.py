#!/usr/bin/python

import sys

import forgi.threedee.model.coarse_grain as ttmc
import borgy.utilities.vector as cuv
import borgy.graph.graph_pdb as cgg
import borgy.utilities.debug as cud

def print_new_bulge_angles(bg):
    '''
    This should be refactored to use the borgy.builder.stats.AngleStat class.
    '''
    for define in bg.defines.keys():
        if define[0] == 'i' or define[0] == 'm':
            (as1, as2) = bg.get_bulge_angle_stats(define)

            print "angle", as1.pdb_name, as1.dim1, as1.dim2, as1.u, as1.v, as1.t, as1.r1, as1.u1, as1.v1, as1.ang_type, " ".join(map(str, as2.define)), " ".join(as1.seqs)
            print "angle", as2.pdb_name, as2.dim1, as2.dim2, as2.u, as2.v, as2.t, as2.r1, as2.u1, as2.v1, as2.ang_type, " ".join(map(str, as2.define)), " ".join(as2.seqs)

        else:
            continue

def print_stem_stats(bg):
    for d in bg.defines.keys():
        if d[0] == 's':

            #if bg.defines[d][0] == 1:
            #    continue

            ss = bg.get_stem_stats(d)

            print "stem", ss.pdb_name, ss.bp_length, ss.phys_length, ss.twist_angle, " ".join(map(str, ss.define))

def get_loop_stat(bg, d):
    base_pair_length = abs(bg.defines[d][0] - bg.defines[d][1])
    phys_length = cuv.magnitude(bg.coords[d][1] - bg.coords[d][0])

    stem1 = list(bg.edges[d])[0]
    stem1_vec = bg.coords[stem1][1] - bg.coords[stem1][0]
    twist1_vec = bg.twists[stem1][1] - bg.twists[stem1][0]
    bulge_vec = bg.coords[d][1] - bg.coords[d][0]

    (r,u,v) = cgg.get_stem_separation_parameters(stem1_vec, twist1_vec, bulge_vec)
    return (base_pair_length, r, u, v)

def print_loop_stats(bg):
    for d in bg.defines.keys():
        if d[0] == 'h':
            (base_pair_length, r, u, v) = get_loop_stat(bg, d)
            print "loop", bg.name, base_pair_length, r, u, v, " ".join(map(str, bg.defines[d]))

def print_5prime_unpaired(bg):
    for d in bg.defines.keys():
        if d[0] == 'f':
            (base_pair_length, r, u, v) = get_loop_stat(bg, d)
            print "5prime", bg.name, base_pair_length, r, u, v, " ".join(map(str, bg.defines[d]))

def print_3prime_unpaired(bg):
    for d in bg.defines.keys():
        if d[0] == 't':
            (base_pair_length, r, u, v) = get_loop_stat(bg, d)
            print "3prime", bg.name, base_pair_length, r, u, v, " ".join(map(str, bg.defines[d]))

def main():
    if len(sys.argv) < 2:
        print >>sys.stderr, "Usage: ./graph_to_angles.py temp.comp"
        print
        print >>sys.stderr, "Traverse a structure and output the stems that are connected by a bulge"
        run_tests()
        sys.exit(1)

    if sys.argv[1] == '-':
        f = sys.stdin
    else:
        f = open(sys.argv[1], 'r')

    bg = ttmc.CoarseGrainRNA(sys.argv[1])
    print_new_bulge_angles(bg)
    print_stem_stats(bg)
    print_loop_stats(bg)
    print_5prime_unpaired(bg)
    print_3prime_unpaired(bg)

if __name__=="__main__":
    main()
