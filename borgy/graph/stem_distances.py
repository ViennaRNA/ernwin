#!/usr/bin/python

import sys, math
from bulge_to_fork_graph import parse_graph, BulgeGraph

def print_neato(bg):
    stems = bg.get_bulged_stem_names()

    # The different nodes for different types of bulges

    node_defs = dict()
    node_lines = dict()

    # loops
    node_defs[1] = '\t{node [shape=circle,color=green,fontsize=12'
    node_defs[2] = '\t{node [shape=circle,color=blue,fontsize=12'
    node_defs[3] = '\t{node [shape=hexagon,color=red,fontsize=12'
    node_defs[4] = '\t{node [shape=octagon,color=red,fontsize=12'
    node_defs[5] = '\t{node [shape=octagon,color=red,fontsize=12'
    node_defs[6] = '\t{node [shape=octagon,color=red,fontsize=12'

    node_lines = ''
    connection_lines = ''

    print "graph G {"
    #print "\tgraph [overlap=scale];"
    print "\tnode [shape=box];"

    for key2 in bg.defines.keys():
        if bg.vert_to_name[key2][0] == 's':
            node_lines += '\t{node [label=\"%s\\n(%d)\"] %s};\n' % (bg.vert_to_name[key2], abs(int(bg.defines[key2][0]) - int(bg.defines[key2][2])) + 1, bg.vert_to_name[key2])
            continue

        node_lines += node_defs[bg.weights[key2]] 
        node_lines += ',label=\"%s \\n(' % (bg.vert_to_name[key2])
        total_bulge = 0

        for j in range(0, len(bg.defines[key2]), 2):
            if j != 0:
                node_lines += ','

            total_bulge += abs((int(bg.defines[key2][j+1]) - int(bg.defines[key2][j]) + 1))
            node_lines += "%d" % (int(bg.defines[key2][j+1]) - int(bg.defines[key2][j]) + 1)
        j = j / 2

        while j < bg.weights[key2]-1:
            node_lines += ",0"
            j += 1

        width = math.sqrt(2.0 * total_bulge / 10.0) 
        height = width

        if bg.weights[key2] == 2:
            node_lines += ")\",width=%f,heigh=%f] %s};\n" % (width, height, bg.vert_to_name[key2])
        else:
            node_lines += ")\"] %s};\n" % (bg.vert_to_name[key2])

    
    for key1 in bg.edges:
        if bg.vert_to_name[key1][0] == 's':
            for key2 in bg.edges[key1]:
                connection_lines += "\t%s -- %s;\n" % (bg.vert_to_name[key1], bg.vert_to_name[key2])

    print node_lines
    print connection_lines
    print "}"


    #print bg.get_named_define

def main():
    if len(sys.argv) < 2:
        print "Usage: ./graph_to_angles.py struct.graph"
        print
        print "Traverse a structure and output the stems that are connected by a bulge"
        sys.exit(1)

    if sys.argv[1] == '-':
        f = sys.stdin
    else:
        f = open(sys.argv[1], 'r')

    bg = parse_graph(f)
    print_neato(bg)

if __name__=="__main__":
    main()
