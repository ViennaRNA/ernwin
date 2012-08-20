#!/usr/bin/python

import sys
from optparse import OptionParser

def replace_first_zeroes(s):
    '''
    Replace all zeroes at the beginning with spaces.
    '''
    news = ''
    i = 0

    while s[i] == '0' and i < len(s):
        news += ' '
        i += 1

    while i < len(s):
        news += s[i]
        i += 1

    return news

def main():
    parser = OptionParser()

    #parser.add_option('-o', '--options', dest='some_option', default='yo', help="Place holder for a real option", type='str')
    #parser.add_option('-u', '--useless', dest='uselesss', default=False, action='store_true', help='Another useless option')

    (options, args) = parser.parse_args()

    ll = 0

    for line in sys.stdin:
        if len(line) > ll:
            ll = len(line)

        print line.strip()

    s = ''
    s10 = ''
    s100 = ''

    for i in range(ll):
        s += str((i + 1) % 10)
        s10  += str(((i + 1) / 10) % 10)
        s100 += str(((i + 1) / 100) % 10)

    print replace_first_zeroes(s100)
    print replace_first_zeroes(s10)
    print replace_first_zeroes(s)


if __name__ == '__main__':
    main()

