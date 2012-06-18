#!/usr/bin/python

import sys
from bulge_to_fork_graph import parse_graph, BulgeGraph
from numpy import array
from graph_to_com_pymol import print_connection

def parse_solution(filename):
    f = open(filename, 'r')

    counter = 1
    solutions = dict()

    for line in f:
        parts = line.strip().split()
        solutions[counter] = array([float(parts[i]) for i in range(len(parts))])
        counter += 1

    return solutions

def get_residue_point_map(bg):
    res_to_point = dict()
    point_counter = 1

    for key in bg.defines.keys():
        #print "key:", key, "define:", bg.defines[key]

        define = bg.defines[key]

        for d in define:
            if int(d) in res_to_point.keys():
                continue
            else:
                res_to_point[int(d)] = point_counter
                point_counter += 1

    return res_to_point

def output_pymol(bg, sol, res_to_point):
    for key in bg.defines.keys():
        define = bg.defines[key]

        if bg.vert_to_name[key][0] == 's':
            print_connection(sol[res_to_point[int(define[0])]], sol[res_to_point[int(define[2])]], 'blue')
            print_connection(sol[res_to_point[int(define[1])]], sol[res_to_point[int(define[3])]], 'blue')

            print_connection(sol[res_to_point[int(define[0])]], sol[res_to_point[int(define[1])]], 'red', 0.4)
            print_connection(sol[res_to_point[int(define[2])]], sol[res_to_point[int(define[3])]], 'red', 0.4)
        else:

            if int(define[0]) > int(define[1]):
                continue
            else:
                try:
                    point1 = res_to_point[int(define[0])-1]
                    point2 = res_to_point[int(define[1])+1]

                    sol1 = sol[point1]
                    sol2 = sol[point2]

                    print_connection(sol1, sol2, 'green')
                except KeyError as ke:
                    continue

def main():
    if len(sys.argv) < 2:
        print "Usage: ./solution_to_pymol.py temp.bulge dg.csv"
        print
        print "Take a bulge graph and the distance geometry solution and create a pymol visualization"
        sys.exit(1)

    if sys.argv[1] == '-':
        f = sys.stdin
    else:
        f = open(sys.argv[1], 'r')

    bg = parse_graph(f)
    sol = parse_solution(sys.argv[2])
    res_to_point = get_residue_point_map(bg)

    print "from pymol.cgo import *"
    print "from pymol import cmd"
    print "from pymol.vfont import plain"
    print "obj1 = ["

    output_pymol(bg, sol, res_to_point)

    print "]"
    print "cmd.load_cgo(obj1, 'ss')"

if __name__=="__main__":
    main()
