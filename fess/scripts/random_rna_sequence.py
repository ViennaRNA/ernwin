#!/usr/bin/python

import sys
import random
from optparse import OptionParser

def main():
    usage = './random_rna_sequence.py [length]'
    usage += 'Create a random RNA sequence of a particular length.'
    parser = OptionParser()

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    length = int(args[0])

    print "".join([random.choice(['A','C','G','U']) for i in range(length)])

if __name__ == '__main__':
    main()

