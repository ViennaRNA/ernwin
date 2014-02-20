#!/usr/bin/python

import collections as c
import itertools as it
import forgi.threedee.utilities.graph_pdb as ftug
import forgi.threedee.utilities.vector as ftuv
import forgi.threedee.model.coarse_grain as ftmc
import forgi.utilities.debug as fud
import numpy as np
import sys
from optparse import OptionParser

def main():
    usage = """
    python average_atom_positions.py file1.pdb file2.pdb ...
    """
    num_args= 0
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    poss = c.defaultdict(list)

    for i,arg in enumerate(args):
        cg = ftmc.from_pdb(arg)
        fud.pv('i')

        for d in cg.defines.keys():
            origin, basis = ftug.element_coord_system(cg, d)

            for i, r in it.izip(it.count(),
                                cg.define_residue_num_iterator(d)):
                atoms = cg.chain[r]
                for a in atoms:
                    avec = a.get_vector().get_array()
                    atom_pos = ftuv.change_basis(avec - origin, basis, ftuv.standard_basis)
                    identifier = "%s %s %d %s" % (d[0], 
                                                  " ".join(map(str, cg.get_node_dimensions(d))),
                                                  i, a.name)
                    poss[identifier] += [atom_pos]

    print "import collections as co"
    print "avg_atom_poss = dict()"

    for key in poss.keys():
        pos = np.mean(poss[key], axis=0)
        print 'avg_stem_poss["%s"] = [%s] #%d' % (key, ",".join(map(str, pos)), len(poss[key]))

if __name__ == '__main__':
    main()

