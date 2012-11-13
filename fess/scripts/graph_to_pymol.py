#!/usr/bin/python

import sys

import corgy.graph.bulge_graph as cgb

def color_stems(bg):
    stems = bg.get_bulged_stem_names()

    for stem_pair in stems:
        #print "stem_pair:", stem_pair

        stem_pair1 = [int(part) for part in bg.get_named_define(stem_pair[0])]
        stem_pair2 = [int(part) for part in bg.get_named_define(stem_pair[1])]

        #print "stem_pair1:", stem_pair1
        #print "stem_pair2:", stem_pair2

        if abs(stem_pair1[0] - stem_pair1[2]) == 0 or abs(stem_pair2[0] - stem_pair2[2]) == 0:
            continue
        else:
            print "color green, resi %d-%d or resi %d-%d" % (stem_pair1[0], stem_pair1[1], stem_pair1[2], stem_pair1[3])
            print "color green, resi %d-%d or resi %d-%d" % (stem_pair2[0], stem_pair2[1], stem_pair2[2], stem_pair2[3])
            if abs(stem_pair2[1] - stem_pair1[3]) > 1:
                print "color green, resi %d-%d" % (stem_pair2[2]+1, stem_pair1[3]-1)
            if abs(stem_pair1[2] - stem_pair2[0]) > 1:
                print "color green, resi %d-%d" % (stem_pair1[1]+1, stem_pair2[0]-1)
        

def hide_stems(bg):
    for key in bg.defines.keys():
        define = bg.defines[key]

        if key[0] == 's':
            #print >>sys.stderr, "key:", key, "define:", define

            for i in range(int(define[0]), int(define[1])+1):
                print "hide cartoon, resi %d" % (i)
                print "hide sticks, resi %d" % (i)

            for i in range(int(define[2]), int(define[3])+1):
                print "hide sticks, resi %d" % (i)
                print "hide cartoon, resi %d" % (i)
        else:
            #print >>sys.stderr,key, "bg.weights:", bg.weights[key]
            if bg.weights[key] == 2:
                for i in range(0, len(define), 2):
                    for j in range(int(define[i])+1,int(define[i+1])):
                        print "color yellow, resi %d" % (j)
            elif bg.weights[key] == 1:
                for i in range(0, len(define), 2):
                    for j in range(int(define[i])+1,int(define[i+1])):
                        print "color blue, resi %d" % (j)
            else:
                for i in range(0, len(define), 2):
                    for j in range(int(define[i])+1,int(define[i+1])):
                        print "color red, resi %d" % (j)
                    

def main():
    if len(sys.argv) < 3:
        print "Usage: ./graph_to_pymol.py struct.graph coarse_grain.py"
        print
        print "Traverse a structure and color the stems that are connected by a bulge"
        sys.exit(1)

    print "hide all"
    print "show sticks, all"
    print "show cartoon, all"
    print "set cartoon_ring_mode"
    print "set cartoon_tube_radius, .3"

    bg = cgb.BulgeGraph(sys.argv[1])
    print "hide all"
    color_stems(bg)
    hide_stems(bg)
    print "run %s" % (sys.argv[2])
    print "show cartoon, temp"
    print "orient"

if __name__=="__main__":
    main()
