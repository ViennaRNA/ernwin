#!/usr/bin/python

import sys, re
import os.path as op
import collections as c
import numpy as np
from optparse import OptionParser

def visit_dir(rmsds, dirname, names):
    if 'out.txt' not in names:
        return

    filename = op.join(dirname, 'out.txt')
    trial_id = op.split(dirname)[1]

    with open(filename) as f:
        file_size = op.getsize(filename) 
        if file_size < 2048:
            return
        f.seek(-1024, 2)
        line = f.readlines()[-1]

        if line.find("native_energy") == 0:
            m = re.search('\[(.*) (.*)\].*min:[ ]+(.*)[ ]+(.*)', line) 

            # group(1) = pdb id
            # group(2) = length in nucleotides
            # group(3) = energy
            # group(4) = rmsd
            rmsds[(int(m.group(2)), m.group(1))] += [(float(m.group(3)), float(m.group(4)), filename, trial_id)]

def summarize_rmsds(rmsds, compact=False, base_dir='', nth=0):
    keys = rmsds.keys()
    #keys.sort(key=lambda x: min([k[1] for k in rmsds[x]]))
    keys.sort()

    min_rmsds = []
    for key1,key2 in keys:
        key = (key1,key2)
        rmsds[key].sort()
        if compact:
            print base_dir, key2, rmsds[key][nth][3]
        else:
            print "[",key2 + "/" + rmsds[key][nth][3],"|", key1,"]", "[", rmsds[key][nth][1], "::", str(rmsds[key][nth][0]) + "]", " ".join([str(k[1]) for k in rmsds[key][1:5]])

        min_rmsds += [rmsds[key][nth][1]]
    min_rmsds.sort()

    if not compact:
        print "average: %.2f median %.2f" % (np.mean(np.array(min_rmsds)), min_rmsds[len(min_rmsds) / 2])

def main():
    usage = './summarize_gibbs_output.py base_dir'
    parser = OptionParser()

    parser.add_option('-n', '--nth-best', dest='nth_best', default=0, help="Display the n-th best energy structure (0-based)", type='int')
    parser.add_option('-c', '--compact', dest='compact', default=False, action='store_true', help='Display a compact representation of the structures consisting of only the name and number.')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    rmsds = c.defaultdict(list)
    for arg in args:
        op.walk(arg, visit_dir, rmsds)

    summarize_rmsds(rmsds, options.compact, args[0], nth=options.nth_best)


if __name__ == '__main__':
    main()

