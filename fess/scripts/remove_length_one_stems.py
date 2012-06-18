#!/usr/bin/python

import sys

def remove_length_one_stems(line):
    stem_pairs = dict()
    opens = []
    line = list(line)
    #print >>sys.stderr, "line:", line, "len:", len(line)

    for i in range(len(line)):
        if line[i] == '(':
            opens.append(i)
        elif line[i] == ')':
            corr = opens.pop()
            stem_pairs[i] = corr
            stem_pairs[corr] = i

    p = '.'
    pp = '.'

    for i in range(len(line)):
        if line[i] == '.':
            if (p == '(' or p == ')') and pp == '.':
                #print >>sys.stderr, "i:", i
                line[i-1] = '.'
                line[stem_pairs[i-1]] = '.'
        pp = p
        p = line[i]

    if pp == '.' and line[i] == ')':
        line[i] = '.'
        line[stem_pairs[i]] = '.'

    return "".join(line)

def main():
    if len(sys.argv) < 2:
        print >>sys.stderr, "Usage: ./remove_length_one_stems.py dotplot"

    f = open(sys.argv[1], 'r')
    line = f.readlines()[0].strip()

    print remove_length_one_stems(line)

if __name__ == "__main__":
    main()
