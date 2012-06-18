#!/usr/bin/python

from optparse import OptionParser
import sys
from corgy.graph.bulge_graph import BulgeGraph
from corgy.visual.pymol import PymolPrinter

def main():
    if len(sys.argv) < 2:
        print "Usage: ./coordinates_to_pymol temp.coordinates"
        sys.exit(1)

    parser = OptionParser()

    parser.add_option('-c', '--centers', dest='centers', default=0, help='Display the centers of each segment.', type='int') 
    parser.add_option('-f', '--flexibility', dest='flex', default=None, help='Location of the flexibility statistics file.', type='string') 
    parser.add_option('-x', '--text', dest='print_text', default=False, action='store_true', help='Print the names of the segments in the pymol output')

    (options, args) = parser.parse_args()
    
    bgs = []
    for arg in args:
        bgs += [BulgeGraph(arg)]

    pymol_printer = PymolPrinter()
    pymol_printer.print_text = options.print_text

    for i in range(len(bgs)):
        if len(bgs) > 1 and i == 0:
            pymol_printer.override_color = 'green'
        elif len(bgs) > 1 and i > 0:
            pymol_printer.override_color = 'red'

        #pymol_printer.override_color = 'red'
        bg = bgs[i]

        if options.centers == 1:
            pymol_printer.centers_to_pymol(bg)
        elif options.flex != None:
            pymol_printer.flex_to_pymol(bg, options.flex)
        
        pymol_printer.coordinates_to_pymol(bg)

    pymol_printer.output_pymol_file()

if __name__ == '__main__':
    main()
