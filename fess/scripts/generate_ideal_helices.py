#!/usr/bin/python

from __future__ import absolute_import
from __future__ import print_function
import sys, os
import os.path as op
import random as rand
import subprocess as sp
import tempfile

from optparse import OptionParser
from six.moves import range

rosetta_base = '/scr/plastilin/pkerp/apps/rosetta/latest/'
rosetta_db = op.join(rosetta_base, 'rosetta_database')
rosetta_helix = op.join(rosetta_base, 'rosetta_source/bin/rna_helix.linuxgccrelease')

def create_seq_fasta_file(stem_length, use_comp=True):
    comp = dict({'a':'u', 'c':'g', 'g':'c', 'u':'a'})

    seq = [rand.choice(['a','c','g','u']) for i in range(stem_length)]

    if use_comp:
        seq += [comp[s] for s in seq[::-1]]

    f = open('stem_seq.fasta', 'w')
    f.write('>seq\n')
    f.write("".join(seq) + "\n")
    f.close()

    return "".join(seq)

def create_helix_using_fiber(seq, stem_length):
    '''
    Create an ideal helix using the fiber program from the X3DNA package.

    @param seq: The sequence of the helix.
    @param stem_length: The length of the helix in nucleotides.
    '''
    ideal_filename = "fess/stats/stems/ideal_%d_%d_%d_%d.pdb" % (1, stem_length, stem_length+1, 2*stem_length)

    tmpdir = tempfile.mkdtemp()
    filename = os.path.join(tmpdir, 'myfifo')

    #ideal_filename = "ideal_%d.pdb" % (stem_length)
    command = ['fiber']
    command += ['-seq=%s' % (seq)]
    command += ['-rna']
    command += [filename]

    print(" ".join(command))
    sp.call(command)

    command = ['sed', '-i', 's/ B/ A/g', filename]
    print(" ".join(command))
    sp.call(command)

    command = ['fess/scripts/make_rna_rosetta_ready.py', filename]
    sp.call(command, stdout=open(ideal_filename, 'w'), stderr=sys.stderr)

    command = ['sed', '-i', 's/rG/ G/g', ideal_filename]
    sp.call(command)
    command = ['sed', '-i', 's/rC/ C/g', ideal_filename]
    sp.call(command)
    command = ['sed', '-i', 's/rA/ A/g', ideal_filename]
    sp.call(command)
    command = ['sed', '-i', 's/rU/ U/g', ideal_filename]
    sp.call(command)

    os.remove(filename)
    os.rmdir(tmpdir)

    return

def create_helix_using_rosetta():
    command = [rosetta_helix]
    command += ["-database", rosetta_db]
    command += ["-fasta", "stem_seq.fasta"]
    command += ["-out:file:silent", "/dev/null"]

    sp.call(command)

def move_helix_to_stem_stats(seq, stem_length):
    ideal_filename = "fess/stats/stems/ideal_%d_%d_%d_%d.pdb" % (1, stem_length, stem_length+1, 2*stem_length)
    command = ["mv"]
    command += [seq+".pdb", ideal_filename]
    sp.call(command)


def main():
    usage = './generate_ideal_helices.py length'
    parser = OptionParser(usage=usage)

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')
    parser.add_option('-p', '--program', dest='program', default='rosetta', help="Which program should we use for generating the ideal helix. The available options are rosetta and fiber", type='str')

    (options, args) = parser.parse_args()

    if len(args) < 1:
        parser.print_help()
        sys.exit(1)

    stem_length = int(args[0])

    if options.program == "rosetta":
        seq = create_seq_fasta_file(stem_length, use_comp = True)
        create_helix_using_rosetta()
        move_helix_to_stem_stats(seq, stem_length)
    else:
        seq = create_seq_fasta_file(stem_length, use_comp = False)
        create_helix_using_fiber(seq, stem_length)

    print() 

if __name__ == '__main__':
    main()

