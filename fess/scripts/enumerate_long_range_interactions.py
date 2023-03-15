#!/usr/bin/python

from __future__ import absolute_import
from __future__ import print_function
import collections as c
import numpy as np
import sys

from optparse import OptionParser

import forgi.threedee.model.coarse_grain as ftmc
from six.moves import range

def create_heatmap(xvals, yvals, zcalculator, filename=None):
    x = np.array(list(range(0, len(xvals)+1)))
    y = np.array(list(range(0, len(yvals)+1)))

    X,Y = np.meshgrid(x, y)
    Z = zcalculator(X,Y)

    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()

    cmap = plt.cm.Blues
    cmap.set_under([0,0,0,0.])
    p = ax.pcolor(X,Y,Z, cmap=cmap, vmin=0.001)

    #ax.set_frame_on(False)
    ax.set_xticks(x + 0.5, minor=False)
    ax.set_yticks(y + 0.5, minor=False)

    ax.set_xticklabels(xvals)
    ax.set_yticklabels(yvals)

    for i in x:
        l = ax.axvline(i-0.010, color='w', ls='-', lw=1.4)
    for i in y:
        ax.axhline(i, color='w', ls='-', lw=1.4)

    for t in ax.xaxis.get_major_ticks(): 
        t.tick1On = False 
        t.tick2On = False 
    for t in ax.yaxis.get_major_ticks(): 
        t.tick1On = False 
        t.tick2On = False 

    locs, labels = plt.xticks()
    plt.setp(labels, rotation=90)

    ax.autoscale_view('tight')

    #plt.colorbar(p)
    fig.set_size_inches(5.0, 4.0)
    if filename == None:
        plt.show()
    else:
        plt.subplots_adjust(left=0.1, right=.95, top=0.97, bottom=0.15)
        plt.savefig(filename)

def main():
    usage = """
    python enumerate_long_range_interactions.py [file1] [file2]

    Go through each file and categorize the list by the types of interactions present.
    """
    num_args= 1
    parser = OptionParser(usage=usage)

    parser.add_option('-o', '--options', dest='filename', default=None, help="Output the file in a particular location", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    interactions = c.defaultdict(set)
    pdb_interactions = c.defaultdict(set)
    all_pdbs = []

    for arg in args:
        cg = ftmc.CoarseGrainRNA(arg)
        all_pdbs += [cg.name]

        for k in cg.longrange.keys():
            for v in cg.longrange[k]:
                interaction_type = [k[0], v[0]]
                interaction_type.sort()
                interaction_type = "".join(interaction_type)
                interactions[interaction_type].add(cg.name)
                pdb_interactions[cg.name].add(interaction_type)

    sorted_pdbs = sorted(all_pdbs, key=lambda x: len(pdb_interactions[x]))
    out_str = "   "
    for p in sorted_pdbs:
        out_str += " %s" % (p)

    interaction_keys = sorted(list(interactions.keys()), key=lambda x: len(interactions[x]))
    print(out_str)

    for i in interaction_keys:
        out_str = "%s:" % (i)
        for p in sorted_pdbs:
            if p in interactions[i]:
                out_str += " ****"
            else:
                out_str += "     "
        print(out_str)


    # Create a heatmap


    @np.vectorize
    def interaction_exists(x, y):
        if y >= len(interaction_keys) or x >= len(sorted_pdbs):
            return 0.
        if sorted_pdbs[x] in interactions[interaction_keys[y]]:
            return 1.
        else:
            return 0.

    create_heatmap(sorted_pdbs, interaction_keys, interaction_exists, options.filename)
    print(len(interaction_keys))

if __name__ == '__main__':
    main()

