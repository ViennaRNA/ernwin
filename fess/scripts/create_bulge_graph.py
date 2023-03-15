#!/usr/bin/python

from __future__ import absolute_import
from __future__ import print_function
import sys, os

import borgy.utilities.debug as cud
from six.moves import range

def error_exit(message):
    print(message, file=sys.stderr)
    sys.exit(1)

# A wrapper for a simple dictionary addition
# Added so that debugging can be amde easier
def add_bulge(bulges, bulge, context, message):
    #print >>sys.stderr,"Adding bulge", context, bulge, message
    #bulge = (context, bulge)
    bulges[context] = bulges.get(context, []) + [bulge]
    return bulges

def any_difference_of_one(stem, bulge):
    '''
    See if there's any difference of one between the two
    ends of the stem [(a,b),(c,d)] and a bulge (e,f)

    @param stem: A couple of couples (2 x 2-tuple) indicating the start and end
                 nucleotides of the stem in the form ((s1, e1), (s2, e2))
    @param bulge: A couple (2-tuple) indicating the first and last position
                  of the bulge.

    @return: True if there is an overlap between the stem nucleotides and the 
                  bulge nucleotides
             False otherwise
    '''
    for stem_part in stem:
        for part in stem_part:
            for bulge_part in bulge:
                if abs(bulge_part - part) == 0:
                    return True
    return False

def create_bulge_graph(stems, bulges):
    '''
    Find out which stems connect to which bulges

    Stems and bulges which share a nucleotide are considered connected.

    @param stems: A list of tuples of tuples of the form [((s1, e1), (s2, e2))]
                  where s1 and e1 are the nucleotides at one end of the stem
                  and s2 and e2 are the nucleotides at the other.

    @param bulges: A list of tuples of the form [(s, e)] where s and e are the 
                   numbers of the nucleotides at the start and end of the bulge.

    @return: A dictionary, indexed by stem names, containing the indexes of the
             the bulges.
    '''
    #print "stems:", stems
    stem_bulges = dict()
    for i in range(len(stems)):
        stem = stems[i]
        stem_str = "s%d" % (i)
        for j in range(len(bulges)):
            bulge = bulges[j]
            #print >>sys.stderr, "stem:", stem, "bulge:", bulge
            if any_difference_of_one(stem, bulge):
                stem_bulges_set = stem_bulges.get(i, set())
                stem_bulges_set.add(j)
                stem_bulges[i] = stem_bulges_set

    #print >>sys.stderr, "stem_bulges:", stem_bulges
    return stem_bulges

# Find out which stems connect to which stems
def create_stem_graph(stems, bulge_counter):
    '''
    Determine which stems are connected to each other. A stem can be connected to
    another stem when there is an interior loop with an unpaired nucleotide on
    one side. In this case, a bulge will be created on the other side, but it
    will only consist of the two paired bases around where the unpaired base 
    would be if it existed.

    The defines for these bulges will be printed as well as the connection strings
    for the stems they are connected to.

    @param stems: A list of tuples of tuples of the form [((s1, e1), (s2, e2))]
                  where s1 and e1 are the nucleotides at one end of the stem
                  and s2 and e2 are the nucleotides at the other.
    @param bulge_counter: The number of bulges that have been encountered so far.

    @return: A dictionary indexed by the number of a stem, containing a set of the 
             other stems that the index is connected to.
    '''
    #print "stems:", stems
    stem_stems = dict()
    define_text = ""
    connect_text = ""
    for i in range(len(stems)):
        stem = stems[i]
        stem_str = "s%d" % (i)
        for j in range(i+1, len(stems)):
            for k1 in range(2):
                # don't fear the for loop
                for k2 in range(2):
                    for l1 in range(2):
                        for l2 in range(2):
                            s1 = stems[i][k1][l1]
                            s2 = stems[j][k2][l2]
                            if abs(s1 - s2) == 1:
                                stem_stems_set = stem_stems.get(i, set())
                                if j not in stem_stems_set:
                                    define_text +=  "define b%d 1 %d %d\n" % (bulge_counter, min(s1, s2)+1, max(s1, s2)+1)
                                    #define_text += "%d %d\n" % (i, j)
                                    connect_text += "connect s%d b%d\n" % (i, bulge_counter)
                                    connect_text += "connect s%d b%d\n" % (j, bulge_counter)
                                    bulge_counter += 1
                                    stem_stems_set.add(j)
                                stem_stems[i] = stem_stems_set
    #print >>sys.stderr, "stem_stems:", stem_stems
    print(define_text)
    print(connect_text)

    return stem_stems


def print_bulge_graph(graph):
    '''
    Print out the connections in the graph.

    @param graph: A dictionary indexed by stem number containing a set
                  of the bulge numbers that it is connected to.
    '''
    for key in graph.keys():
        stem_str = "connect s%d" % (key)
        for item in graph[key]:
            stem_str += " b%d" % (item)
        print(stem_str)

def print_stems(stems):
    '''
    Print the names and definitions of the stems.

    @param stems: A list of tuples of tuples of the form [((s1, e1), (s2, e2))]
                  where s1 and e1 are the nucleotides at one end of the stem
                  and s2 and e2 are the nucleotides at the other.
    '''
    for i in range(len(stems)):
        # one is added to each coordinate to make up for the fact that residues are 1-based
        ss1 = stems[i][0][0]+1
        ss2 = stems[i][0][1]+1
        se1 = stems[i][1][0]+1
        se2 = stems[i][1][1]+1

        print("define s%d 0 %d %d %d %d" % (i, min(ss1,se1), max(ss1,se1), min(ss2,se2), max(ss2,se2)))

def print_bulges(bulges):
    '''
    Print the names and definitions of the bulges.

    @param bulges: A list of tuples of the form [(s, e)] where s and e are the 
                   numbers of the nucleotides at the start and end of the bulge.
    '''
    for i in range(len(bulges)):
            #print "bulge:", bulge
        bulge_str = "define b%d 1" % (i)
        bulge = bulges[i]
        #print >>sys.stderr, "bulge:", bulge
        bulge_str += " %d %d" % (bulge[0]+1, bulge[1]+1)
        print(bulge_str)

def condense_stem_pairs(stem_pairs):
    '''
    Given a list of stem pairs, condense them into stem definitions

    I.e. the pairs (0,10),(1,9),(2,8),(3,7) can be condensed into
    just the ends of the stem: [(0,10),(3,7)]

    @param stem_pairs: A list of tuples containing paired base numbers.

    @return: A list of tuples of tuples of the form [((s1, e1), (s2, e2))]
                  where s1 and e1 are the nucleotides at one end of the stem
                  and s2 and e2 are the nucleotides at the other.
    '''
    stem_pairs.sort()

    prev_pair = (-10, -10)

    stems = []
    start_pair = None

    for pair in stem_pairs:
        # There's a potential bug here since we don't check the direction
        # but hopefully it won't bite us in the ass later
        if abs(pair[0] - prev_pair[0]) != 1 or abs(pair[1] - prev_pair[1]) != 1:
            if start_pair != None:
                stems += [(start_pair, prev_pair)]
            start_pair = pair
    
        prev_pair = pair

    if start_pair != None:
        stems += [(start_pair, prev_pair)]

    return stems

def print_brackets(brackets):
    '''
    Print the brackets and a numbering, for debugging purposes

    @param brackets: A string with the dotplot passed as input to this script.
    '''
    numbers = [chr(ord('0') + i % 10) for i in range(len(brackets))]
    tens = [chr(ord('0') + i / 10) for i in range(len(brackets))]
    print("brackets:\n", brackets, "\n", "".join(tens), "\n" ,"".join(numbers))

def find_bulges_and_stems(brackets):
    '''
    Iterate through the structure and enumerate the bulges and the stems that are
    present.

    The returned stems are of the form [[(s1, s2), (e1,e2)], [(s1,s2),(e1,e2)],...]
    where (s1,s2) are the residue numbers of one end of the stem and (e1,e2) are the
    residue numbers at the other end of the stem
    (see condense_stem_pairs)

    The returned bulges are of the form [(s,e), (s,e),...] where s is the start of a bulge
    and e is the end of a bulge

    @param brackets: A string with the dotplot passed as input to this script.
    '''
    prev = '.'
    context = 0

    bulges = dict()
    finished_bulges = []
    context_depths = dict()

    opens = []
    stem_pairs = []

    stems = dict()

    dots_start = 0
    dots_end = 0

    context_depths[0] = 0
    cud.pv('len([b for b in brackets if b == "("])')
    cud.pv('len([b for b in brackets if b == ")"])')

    cud.pv('brackets')

    i = 0
    for i in range(len(brackets)):
        #print >> sys.stderr, "bracket:", brackets[i]
        if brackets[i] == '(':
            opens.append(i)

            #print >>sys.stderr,i+1, "Open bracket, context:", context, "depth", context_depths[context]

            if prev == '(':
                context_depths[context] = context_depths.get(context, 0) + 1
                continue
            else:
                context += 1
                context_depths[context] = 1

            if prev == '.':
                dots_end = i-1
                bulges = add_bulge(bulges, (dots_start-1, dots_end+1), context, "4")
            """
            if prev == ')':
                bulges = add_bulge(bulges, (-1,-1), context, "3")
            """

        if brackets[i] == ')':
            if len(opens) == 0:
                error_exit("ERROR: Unmatched close bracket")

            stem_pairs.append((opens.pop(), i))
            #print >>sys.stderr,i+1, "Close bracket, context:", context, "depth:", context_depths[context]

            if prev == '(':
                bulges = add_bulge(bulges, (i-1, i), context, "4")

            context_depths[context] -= 1

            if context_depths[context] == 0:
                cud.pv('bulges')
                finished_bulges += bulges[context]
                bulges[context] = []
                context -= 1
 

            if prev == '.':
                dots_end = i-1
                bulges = add_bulge(bulges, (dots_start-1, dots_end+1), context, "2")
  
        if brackets[i] == '.':
            if prev == '.':
                continue

            dots_start = i

        prev = brackets[i]
        #print i, "bulges:",context, bulges.get(context, [])
    if prev == '.':
        dots_end = i
        bulges = add_bulge(bulges, (dots_start-1, dots_end), context, "7")
    elif prev == '(':
        print("Unmatched bracket at the end", file=sys.stderr)
        sys.exit(1)
    """
    elif prev == ')':
        bulges = add_bulge(bulges, (i+1, i), context, "8")
    """
    
    if context in bulges.keys():
        finished_bulges += bulges[context]

    if len(opens) > 0:
        error_exit("ERROR: Unmatched open bracket")

    #print >>sys.stderr, "finished_bulges:", finished_bulges

    stem_pairs.sort()
    #print "End, finished_bulges:", finished_bulges
    #print "stem_pairs:", stem_pairs
    stems = condense_stem_pairs(stem_pairs)
    #print "stems:", stems
    
    return (finished_bulges, stems)

def print_name(filename):
    print("name", os.path.splitext(filename)[0])

def main():
    if len(sys.argv) < 2:
        print("""
        Usage: ./create_bulge_graph.py file.dotbracket [name]"
        
        Creates a graph of the paired and unpaired regions within a
        dotplot. Paired regions are called stems while unpaired regions
        are called bulges.

        The file created contains four sections:
        1. Name (optional):
            The name of this graph. This can be used, for example, to specify which
            PDB file the dotplot was inferred from. The name should not contain any
            spaces.

        2. Length (optional):
            The length of the dotplot that was used to create this graph.

        3. Definitions:
            Each line in the definitions sections begins with the keyword 'define'
            It is followed by the name of the region and the numbers of the 
            nucleotides which define that region.

            The numbers of the nucleotides always come in sets of two for bulges and
            sets of four for stems (which can be interpreted as two sets of two). Each
            set of two numbers indicates the starting and ending nucleotide for that
            region.

            All numbers are 1-based.

        4. Connections:
            This section shows how the different regions are connected to each other.
            The term connected means that they share a common nucleotide.

            Each line begins with the keyword 'connect'. It is followed by the name
            of the region for which the connections are being described. The name
            is followed by a number of different names. Each of these regions is connected
            to the first.

            One region may be defined by more than one connect statement. This section
            of the file is simpy an adjacency list.


        The results are printed to standard out.
        """)
        sys.exit(1)
    if sys.argv[1] == '-':
        f = sys.stdin

    f = open(sys.argv[1])
    brackets = "".join(f.readlines()).replace('\n', '')
        #print "stems:", stems

    if len(sys.argv) == 3:
        print("name", sys.argv[2])
        print("length", len(brackets))

    (bulges, stems) = find_bulges_and_stems(brackets)

    print_stems(stems)
    print_bulges(bulges)

    bulge_graph = create_bulge_graph(stems, bulges)
    stem_graph = create_stem_graph(stems, len(bulges))

    print_bulge_graph(bulge_graph)

if __name__ == "__main__":
    main()
