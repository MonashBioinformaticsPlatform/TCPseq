#!/usr/bin/env python
import getopt
import sys
import gzip

options, remainder = getopt.getopt(sys.argv[1:], 'i:f:',
                                   ['input=',
                                   'fastq='])

for opt, arg in options:
    if opt in ('-i', '--input'):
        input_sam = arg
    if opt in ('-f', '--fastq'):
        input_fq = arg
        

to_rescue=dict()
counter = 0
millions = 0
with open(input_sam, 'r') as fin:
    curr_read = "x"
    forward = False
    mapped = False
    primary = ""
    for line in fin:
        splut = line.rstrip('\n').split('\t')
        if splut[0].startswith("@"):
            continue
        counter += 1
        if counter % 1000000 == 0:
            millions += 1
            sys.stderr.write("%d million alignments done\n" % (millions))
        if splut[0] != curr_read:
            if not mapped:
                # print the read to STDOUT if it never mapped at all
                if not primary == "":
                    print primary
            elif not forward:
                # flag the read to extract from the original fastq file, if it only mapped in a reverse orientation
                # note: original pipeline took this from the SAM file, but this is entry will be the rev-comped alignment, not the read itself
                to_rescue['@'+curr_read] = 'y'
            forward = False
            mapped = False
            primary = ""
            curr_read = splut[0]
        flag = int(splut[1])
        if not (flag & 0x10):  # if not (rev comp)
            forward = True
        if not (flag & 0x4):  # if not (unmapped)
            mapped = True
        if not (flag & 0x100):  # if not (secondary alignment)
            primary = "\n".join([ "@"+splut[0], splut[9], "+", splut[10] ])
# print the very last read if it never mapped in a forward orientation
if not mapped:
    print primary
elif not forward:
   to_rescue[curr_read] = 'y'
   
sys.stderr.write( "length of to_rescue: \n")
sys.stderr.write( "%d \n" % (len(to_rescue)) )

counter=0
for x in to_rescue:
    sys.stderr.write("%s \n" % (x))
    counter+=1
    if counter > 10:
        break

counter=0
countdown=0
millions = 0
with gzip.open(input_fq, 'r') as fin:
    for line in fin:
        counter += 1
        if counter % 4 == 1:
            rname = line.rstrip().split(' ')
            if counter % 4000000 == 1:
                sys.stderr.write("%s\n" % (rname))
            if rname[0] in to_rescue:
                countdown = 4
        if countdown > 0:
            countdown -= 1
            print line.rstrip()
        if counter % 1000000 == 0:
            millions += 1
            sys.stderr.write("%d million reads done\n" % (millions))
