#!/usr/bin/python

import sys

from borgy.utilities.mcannotate import iterate_over_interactions, iterate_over_residue_list, parse_chain_base, get_dotplot


def main():
    """docstring for main"""
    if len(sys.argv) < 2:
        print "Usage: ./mcannotate_to_dotplot.py"
        exit(1)

    lines = open(sys.argv[1], 'r').readlines()
    
    print get_dotplot(lines)

    pass

if __name__ == '__main__':
    main()
