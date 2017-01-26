#!/usr/bin/env python
import getopt
import sys
options, remainder = getopt.getopt(sys.argv[1:], 'i:', ['input='])

for opt, arg in options:
    if opt in ('-i', '--input'):
        input_sam = arg

with open(input_sam, 'r') as fin:
    curr_read = "x"
    forward = False
    primary = ""
    for line in fin:
        splut = line.rstrip('\n').split('\t')
        if splut[0].startswith("@"):
            continue
        if splut[0] != curr_read:
            if not forward and primary != "":
                # print the read to STDOUT if it never mapped in a forward orientation
                print primary
            forward = False
            primary = ""
            curr_read = splut[0]
        flag = int(splut[1])
        if not ((flag & 0x4) | (flag & 0x10)):  # if not (unmapped or rev conp)
            forward = True
        if not (flag & 0x100):  # if not (secondary alignment)
            primary = "\n".join([ "@"+splut[0], splut[9], "+", splut[10] ])
# print the very last read if it never mapped in a forward orientation
if not forward and primary != "":
    print primary
