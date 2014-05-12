#!/usr/bin/python

import os.path as op
import sys

import forgi.threedee.model.coarse_grain as ttmc
import forgi.threedee.utilities.vector as ftuv
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.utilities.debug as fud

from optparse import OptionParser

def print_new_bulge_angles(bg):
    '''
    This should be refactored to use the borgy.builder.stats.AngleStat class.
    '''
    out_str = ''
    for define in bg.defines.keys():
        if define[0] == 'i' or define[0] == 'm':
            (as1, as2) = bg.get_bulge_angle_stats(define)

            out_str += "angle " + " ".join(map(str, [ as1.pdb_name, as1.dim1, as1.dim2, as1.u, as1.v, as1.t, as1.r1, as1.u1, as1.v1, as1.ang_type, " ".join(map(str, as2.define)), " ".join(as1.seqs), "\n"]))
            out_str += "angle " + " ".join(map(str, [as2.pdb_name, as2.dim1, as2.dim2, as2.u, as2.v, as2.t, as2.r1, as2.u1, as2.v1, as2.ang_type, " ".join(map(str, as2.define)), " ".join(as2.seqs), "\n"]))

        else:
            continue

    return out_str

def print_stem_stats(bg):
    out_str = ''
    for d in bg.defines.keys():
        if d[0] == 's':

            #if bg.defines[d][0] == 1:
            #    continue

            ss = bg.get_stem_stats(d)

            out_str += "stem " + " ".join(map(str, [ ss.pdb_name, ss.bp_length, ss.phys_length, ss.twist_angle, " ".join(map(str, ss.define)), "\n"]))

    return out_str

def get_loop_stat(bg, d):
    base_pair_length = bg.get_length(d) #abs(bg.defines[d][0] - bg.defines[d][1])
    phys_length = ftuv.magnitude(bg.coords[d][1] - bg.coords[d][0])

    stem1 = list(bg.edges[d])[0]
    (s1b, s1e) = bg.get_sides(stem1, d)

    stem1_vec = bg.coords[stem1][s1b] - bg.coords[stem1][s1e]
    twist1_vec = bg.twists[stem1][s1b] - bg.twists[stem1][s1e]
    bulge_vec = bg.coords[d][1] - bg.coords[d][0]

    (r,u,v) = ftug.get_stem_separation_parameters(stem1_vec, twist1_vec, bulge_vec)
    return (base_pair_length, r, u, v)

def print_loop_stats(bg):
    out_str = ''
    for d in bg.defines.keys():
        if d[0] == 'h':
            (base_pair_length, r, u, v) = get_loop_stat(bg, d)
            out_str += "loop " + " ".join(map(str, [ bg.name, base_pair_length, r, u, v, " ".join(map(str, bg.defines[d])), "\n"]))

    return out_str

def print_5prime_unpaired(bg):
    out_str = ''
    for d in bg.defines.keys():
        if d[0] == 'f':
            (base_pair_length, r, u, v) = get_loop_stat(bg, d)
            out_str += "5prime " + " ".join(map(str, [ bg.name, base_pair_length, r, u, v, " ".join(map(str, bg.defines[d])), "\n"]))

    return out_str

def print_3prime_unpaired(bg):
    out_str = ''
    for d in bg.defines.keys():
        if d[0] == 't':
            (base_pair_length, r, u, v) = get_loop_stat(bg, d)
            out_str += "3prime " + " ".join(map(str, [bg.name, base_pair_length, r, u, v, " ".join(map(str, bg.defines[d])), "\n"]))

    return out_str

def main():
    if len(sys.argv) < 2:
        print >>sys.stderr, "Usage: ./graph_to_angles.py temp.comp"
        print
        print >>sys.stderr, "Traverse a structure and output the stems that are connected by a bulge"
        run_tests()
        sys.exit(1)

    parser = OptionParser()

    parser.add_option('-d', '--dump', dest='dump', 
                      default=False, help="Dump a file called temp.angles in the same directory as the original file",
                      action="store_true")

    (options, args) = parser.parse_args()

    bg = ttmc.CoarseGrainRNA(args[0])
    out_str = print_new_bulge_angles(bg)
    out_str += print_stem_stats(bg)
    out_str += print_loop_stats(bg)
    out_str += print_5prime_unpaired(bg)
    out_str += print_3prime_unpaired(bg)

    if options.dump:
        # dump an output file
        if args[0] == '-':
            dirname = '.'
        else:
            dirname = op.dirname(args[0])

        with open(op.join(dirname, 'temp.angles'), 'w') as f:
            f.write(out_str)
    else:
        print out_str

if __name__=="__main__":
    main()
