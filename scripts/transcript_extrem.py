#!/usr/bin/env python

###############################
# ESTIMATES LIBERAL UTR LIMITS BASED ON NEXT ORF START / END

import getopt
import sys
import os
import re
from Bio import SeqIO

NGL_5p = "annotation/Nagal_5pUTRs.bed"
NGL_3p = "annotation/Nagal_3pUTRs.bed"
GBK_DIR = "genbank"
FPUTR_INTRONS = "annotation/fp_UTR_introns.txt"
MANUALLY_FLAGGED = list()

options, remainder = getopt.getopt(sys.argv[1:], 'f:t:g:i:b:', ['fivep=', 
                                                         'threep=',
                                                         'gbk=',
                                                         'introns_5p=',
                                                         'blacklist='
                                                         ])
for opt, arg in options:
    if opt in ('-f', '--fivep'):
        NGL_5p = arg
    elif opt in ('-t', '--threep'):
        NGL_3p = arg
    elif opt in ('-g', '--gbk'):
        GBK_DIR = arg
    elif opt in ('-i', '--introns_5p'):
         FPUTR_INTRONS = arg
    elif opt in ('-b', '--blacklist'):
         MANUALLY_FLAGGED = arg.split(" ")      
gb_files = os.listdir(GBK_DIR)

lastpos = 0
overlap = int()
nonoverlap = int()
HoA = dict()

for gbfile in gb_files:
    for seq_record in SeqIO.parse(GBK_DIR+'/'+gbfile, "genbank"):
        CDSs = []
        chr = "chr" + re.sub("EF4", "", seq_record.id)
        for feature in seq_record.features:
            if feature.type == "CDS":
                CDSs.append(feature)
                if feature.location.start < lastpos:
                    overlap += 1
                else:
                    nonoverlap += 1
                lastpos = feature.location.end
        
        #assumes already ordered by either 3' or 5' ends
        for i in range(0,len(CDSs)):
            gene = CDSs[i].qualifiers['gene'][0]
            line = [gene, int(CDSs[i].location.start), int(CDSs[i].location.end), CDSs[i].location.strand, chr]
            # find upstream gene end: work from i-1 to 0
            up_lim = 0
            if i > 0:
                for j in range(1,i):
                    if CDSs[i-j].location.end < CDSs[i].location.start:
                        up_lim = int(CDSs[i-j].location.end)
                        break   
            #  Find downstream gene start:
            down_lim = len(seq_record)
            if i < len(CDSs)-1:
                for j in range(i+1, len(CDSs)-1):
                    if CDSs[j].location.start > CDSs[i].location.end:
                        down_lim = int(CDSs[j].location.start)
                        break     
            line.append(up_lim)
            line.append(down_lim)
            cline = "\t".join(map(str, line))
            if int(line[3]) == 1:
                utr_5p = line[1] - line[5]
                utr_3p = line[6] - line[2]
            else:
                utr_3p = line[1] - line[5]
                utr_5p = line[6] - line[2]
            line.append(utr_5p)
            line.append(utr_3p)
            line.append(line[2]-line[1])
            if gene in MANUALLY_FLAGGED:
                manflag = 'TRUE'
            else:
                manflag = 'FALSE'
            HoA[gene]=line
            HoA[gene].extend(['NA','NA', 'FALSE', manflag])

            #   1           2           3           4           5           6         
            #   CDS-start  CDS-end      strand      chr         up_lim      down_lim   
#print "overlap:  ", overlap
#print "non-overlapping:  ", nonoverlap

## ALSO INCLUDE NAGALAKSHMI DATA
with open(NGL_5p, 'r') as fin5p:
    for line in fin5p:
        splut = line.rstrip('\n').split('\t')
        if len(splut) == 6:
            genename = splut[3].replace('_5UTR', '')
            if genename in HoA:
                HoA[genename][10] = str(int(splut[2]) - int(splut[1]))
with open(NGL_3p, 'r') as fin3p:
    for line in fin3p:
        splut = line.rstrip('\n').split('\t')
        if len(splut) == 6:
            genename = splut[3].replace('_3UTR', '')
            if genename in HoA:
                HoA[genename][11] = str(int(splut[2]) - int(splut[1]))
                
## flag 5'UTR-intron genes (Mostly Zhang et al 2008; produced from batch download from http://yeastmine.yeastgenome.org querying 'FivePrimeUTRIntron'; 17/8/2016)
with open(FPUTR_INTRONS, 'r') as fin_introns5p:
    for line in fin_introns5p:
        gene = line.rstrip()
        if gene in HoA:
            HoA[gene][12] = 'TRUE'

## Print to STDOUT
print "\t".join(['gene', 'start','end','strand','chr', 'upstr_adj', 'dnstr_adj', 'max_5p', 'max_3p', 'orf_len', 'real_5p', 'real_3p', 'intron_5pUTR', 'manually_flagged'])
for g in HoA:
    print "\t".join([str(x) for x in HoA[g]])
