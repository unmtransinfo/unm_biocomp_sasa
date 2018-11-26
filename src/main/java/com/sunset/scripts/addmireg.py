#!/usr/bin/env python

import sys

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print "usage: python addmireg.py input.rdf output.rdf start"
        sys.exit(1)
    input = open(sys.argv[1],'r')
    output = open(sys.argv[2],'w')
    start = int(sys.argv[3])
    for line in input:
        if line.startswith('$MFMT'):
            output.write('$MFMT $MIREG ' + str(start) + '\n')
            start = start + 1
        else:
            output.write(line)
    output.close()
    input.close()