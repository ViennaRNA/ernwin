#!/usr/bin/python

from __future__ import absolute_import
from __future__ import print_function
import json
import sys
from optparse import OptionParser
import logging
import forgi.utilities.debug as fud
import forgi.graph.bulge_graph as fgb
from six.moves import map

log = logging.getLogger(__name__)
class MotifEntry:
    '''
    An entry in a motif.

    This records the pdb file that this motif came from as well
    as the nucleotides it contains.
    '''

    def __init__(self, line):
        # the define should indicate starting and ending nucleotides
        # of contiguous regions in a alternating fashion:
        # start1 end1 start2 end2
        self.define = []
        self.parse_line(line)

    def parse_line(self, line):
        '''
        Parse each line of the motif atlas.

        Example of a line:

        "2ZI0|1|D|U|9","2ZI0|1|D|U|10","2ZI0|1|D|A|11","2ZI0|1|C|U|9","2ZI0|1|C|U|10","2ZI0|1|C|A|11"

        '''
        parts = line.split(",")
        pdb_id = parts[0].split('|')[0]
        print(("MotifEntry", parts))

        prev_resnum = -10

        for part in parts:
            sparts = part.split('|')

            # find contiguous regions of nucleotides
            resnum = fgb.RESID("{}:{}".format(parts[2],int(parts[4].strip('"'))))
            self.define += [resnum]

        print((line, self.define))

class MotifAlignment:
    '''
    An alignment within a motif.
    '''
    def __init__(self, json_alignment, chainbreak):
        self.residues = []
        self.json_alignment = json_alignment
        self.chainbreak = int(chainbreak)
        self.struct = ''
        self.chains = set()

        self.parse_lines(json_alignment)

    def parse_lines(self, json_alignment):
        '''
        Parse the lines of a motif alignment.

        Example lines:

            ['"3V2F|1|A|U|2687","3V2F|1|A|U|2688","3V2F|1|A|C|2690",
            "3V2F|1|A|C|2691","3V2F|1|A|G|2718","3V2F|1|A|G|2719",
            "3V2F|1|A|U|2720","3V2F|1|A|A|2721","3V2F|1|A|G|2722"']

        @param lines: The lines in the csv representation of the motif atlas.
        '''
        for part in json_alignment:
            sparts = part.split('|')
            self.struct = sparts[0]
            self.model = sparts[1]
            self.chains.add(sparts[2])
            #print "sparts:", sparts

            resnum = fgb.RESID("{}:{}".format(sparts[2],int(sparts[4].strip('"'))))
            self.residues += [resnum]

    def __str__(self):
        return "%s %s | %s" % (self.struct,
                               " ".join(map(str, self.residues[:self.chainbreak])),
                               " ".join(map(str, self.residues[self.chainbreak:])))

class Motif:
    '''
    A motif from the BGSU motif atlas.
    '''

    def __init__(self, json_motif):
        '''
        Initialize a motif using the lines from the csv representation
        of the motif atlas
        '''
        self.entries = []
        self.json_motif = json_motif
        self.alignments = []
        self.chainbreak = None

        self.parse_alignments(self.json_motif)

        pass

    def parse_alignment(self, json_motif):
        '''
        Go through all the alignments in this motif and parse
        them for their source structure as well as involved nucleotides.

        @param json_motif: The section of the motif atlas json file
                           corresponding to this motif.
        @return: Nothing, just fill in the entries in this data structure.
        '''
        if json_motif['alignment']:
            for a in motif['alignment']:
                self.alignments += MotifAlignment(a, json_motif['chainbreak'])

class MotifAtlas:
    '''
    A class to hold a representation of the BGSU RNA 3D Motif Atlas.
    '''

    def __init__(self, filename):
        '''
        Create the class. All it needs is the filename containing
        the motif atlas.

        @param filename: The csv file containing the motif atlas.
        '''
        self.filename = filename
        self.motifs = dict()

        self.load_motif_atlas(filename)

    def load_motif_atlas(self, filename):
        # Load the motif atlas by parsing the csv file
        prev_lines = []
        motif_id = ''

        with open(filename, 'r') as f:
            s = "".join(f.readlines())
            ma = json.loads(s)

            if len(ma) == 0:
                log.warning("Loaded empty motif atlas... something probably went wrong.")

        for motif in ma:
            motif_id = motif["motif_id"].split('.')[0]
            self.motifs[motif_id] = motif

    def struct_in_motif(self, motif_id, struct_name):
        '''
        Check if a particular structure is in themotif_entry alignment
        of the motif.

        @param motif_id: The id of the motif.
        @param struct_name: The name of the structure.
        @return: True or False
        '''
        struct_name = struct_name.upper()

        motif = self.motifs[motif_id]
        if motif['alignment']:
            for a in motif['alignment']:
                if a.find(struct_name) >= 0:
                    return True

        return False

def main():
    usage = """
    python motif_atlas.py internal_loop_motif_atlas.json
    """
    num_args= 1
    parser = OptionParser(usage=usage)

    parser.add_option('-s', '--struct', dest='struct', default=None, help="List the motifs that contain a particular structure", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    if len(args) < num_args:
        parser.print_help()
        sys.exit(1)

    ma = MotifAtlas(args[0])
    if options.struct is not None:
        # search for a particular structure
        for motif in ma.motifs.values():
            if motif['alignment']:
                for a in motif['alignment']:
                    if a.find(options.struct.upper()) >= 0:
                        print((motif["motif_id"]))
                        al = MotifAlignment(motif['alignment'][a], motif['chainbreak'])
                        print (al)

if __name__ == '__main__':
    main()
