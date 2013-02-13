#!/usr/bin/python

import sys

import corgy.graph.bulge_graph as cgb
import corgy.utilities.vector as cuv
import corgy.graph.graph_pdb as cgg

def print_new_bulge_angles(bg):
    '''
    This should be refactored to use the corgy.builder.stats.AngleStat class.
    '''
    for define in bg.defines.keys():
        if define[0] != 's' and len(bg.edges[define]) == 2:
            (as1, as2) = bg.get_bulge_angle_stats(define)

            print "angle", as1.pdb_name, as1.dim1, as1.dim2, as1.u, as1.v, as1.t, as1.r1, as1.u1, as1.v1, as1.ang_type, " ".join(map(str, as2.define))
            print "angle", as2.pdb_name, as2.dim1, as2.dim2, as2.u, as2.v, as2.t, as2.r1, as2.u1, as2.v1, as2.ang_type, " ".join(map(str, as2.define))

        else:
            continue

def print_stem_stats(bg):
    for d in bg.defines.keys():
        if d[0] == 's':
            ss = bg.get_stem_stats(d)

            print "stem", ss.pdb_name, ss.bp_length, ss.phys_length, ss.twist_angle, " ".join(map(str, ss.define))

def print_loop_stats(bg):
    for d in bg.defines.keys():
        if d[0] != 's':
            if bg.weights[d] == 1 and len(bg.edges[d]) == 1 and bg.defines[d][0] != 0 and bg.defines[d][1] != bg.length:
                # unpaired region
                base_pair_length = abs(bg.defines[d][0] - bg.defines[d][1])
                phys_length = cuv.magnitude(bg.coords[d][1] - bg.coords[d][0])
                
                print "loop", bg.name, base_pair_length, phys_length, " ".join(map(str, bg.defines[d]))

def print_5prime_unpaired(bg):
    for d in bg.defines.keys():
        if d[0] != 's':
            if len(bg.edges[d]) == 1 and bg.defines[d][0] == 0:
                base_pair_length = abs(bg.defines[d][0] - bg.defines[d][1])
                phys_length = cuv.magnitude(bg.coords[d][1] - bg.coords[d][0])
                
                print "5prime", bg.name, base_pair_length, phys_length, " ".join(map(str, bg.defines[d]))

def print_3prime_unpaired(bg):
    for d in bg.defines.keys():
        if d[0] != 's':
            if len(bg.edges[d]) == 1 and bg.defines[d][1] == bg.length:
                base_pair_length = abs(bg.defines[d][0] - bg.defines[d][1])
                phys_length = cuv.magnitude(bg.coords[d][1] - bg.coords[d][0])
                
                print "3prime", bg.name, base_pair_length, phys_length, " ".join(map(str, bg.defines[d]))

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

    bg = cgb.BulgeGraph(sys.argv[1])
    print_new_bulge_angles(bg)
    print_stem_stats(bg)
    print_loop_stats(bg)
    print_5prime_unpaired(bg)
    print_3prime_unpaired(bg)

if __name__=="__main__":
    main()
