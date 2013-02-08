#!/usr/bin/python

import sys
import re
import collections as c
from optparse import OptionParser

def main():
    usage = 'find best | grep out.txt | xargs -n 1 tail -n 1 | ./summarize_gibbs_output.py'
    parser = OptionParser()

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < 0:
        parser.print_help()
        sys.exit(1)

    rmsds = c.defaultdict(list)

    for line in sys.stdin:
        if line.find("native_energy") == 0:
            m = re.search('\[(.*) (.*)\].*min:[ ]+(.*)[ ]+(.*)', line) 
            '''
            print m.groups()
            parts = line.strip().split(' ')
            print parts
            '''
            rmsds[(int(m.group(2)), m.group(1))] += [(float(m.group(3)), float(m.group(4)))]

    keys = rmsds.keys()
    #keys.sort(key=lambda x: min([k[1] for k in rmsds[x]]))
    keys.sort()

    for key1,key2 in keys:
        key = (key1,key2)
        rmsds[key].sort()
        print "[",key2, key1,"]", "[", rmsds[key][0][1], "]", " ".join([str(k[1]) for k in rmsds[key][1:5]])

if __name__ == '__main__':
    main()

