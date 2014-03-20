#!/usr/bin/python

import sys
import random
from optparse import OptionParser

def main():
    usage = './random_rna_sequence.py [length]'
    usage += 'Create a random RNA sequence of a particular length.'
    parser = OptionParser()

    parser.add_option('-o', '--output-file', dest='output_file', default=None, help="The filename to output the string to", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    length = int(args[0])

    if options.output_file == None:
        output_file = sys.stdout
    else:
        output_file = open(options.output_file, 'w')

    output_file.write("".join([random.choice(['A','C','G','U']) for i in range(length)]))
    output_file.write("\n")
    output_file.flush()
    output_file.close()

if __name__ == '__main__':
    main()

