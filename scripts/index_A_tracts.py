#!/usr/bin/env python

# adds a flag to table output of sam2table_v2.pl
# e.g python index_A_tracts.py

import re, glob, sys, getopt, gzip

flank_len = 1000
mRNAs_fasta = '../sacCer3/mRNAs_1000flank/mRNAs_1000flank.fa'
#input_filenames = 'tabulated_data/noRTS_v_orf1000exdub_*_m17al.txt.gz'
tbl_fns_pref = 'tabulated_data/tbl1_'
tbl_fns_suff = '_min17.txt.gz'
output_pref = 'tbl2_'
output_suff = '.bed.gz'

options, remainder = getopt.getopt(sys.argv[1:], 'p:h:s:m:f:o:', ['prefix=',
                                                                'handles=',
                                                                'suffix=',
                                                                'mRNA=',
                                                                'flanking=',
                                                                'output_prefix'
                                                                ])
for opt, arg in options:
    if opt in ('-p', '--prefix'):
        tbl_fns_pref = arg
    elif opt in ('-h', '--handles'):
        tbl_fns_uniq = arg.split(':::')
    elif opt in ('-s', '--suffix'):
        tbl_fns_suff = arg
    elif opt in ('-m', '--mRNA'):
        mRNAs_fasta = arg
    elif opt in ('-f', '--flanking' ):
        flank_len = int(arg)
    elif opt in ('-o', '--output_prefix'):
        output_pref = arg

#tbl_fns = [tbl_fns_pref+x+tbl_fns_suff for x in tbl_fns_uniq]
tbl_fns = dict()  #key: in file value: out file
for x in tbl_fns_uniq:
    tbl_fns[tbl_fns_pref+x+tbl_fns_suff] = output_pref+x+output_suff

###########################################
# Store all pA-tract positions in dict of lists 
###########################################

#with open('/Users/sarcher/VMs/bioinf/shared_data/common_data/genomes/sacCer3/mR_900flank.fa', 'r') as fin:
with open(mRNAs_fasta, 'r') as fin:   
    gene = 'First'
    seq = ''
    gene_num = 0
    wg_index = dict()
    gene_len = dict()
    for line in fin:
        line = line.rstrip('\n')
        mobj = re.match('^>\s*(.+)', line)
        if mobj:
            if not gene == "First":
                A_count = 0
                A_tracts = [0] * len(seq)
                for i in range(len(seq)-1, -1, -1): # work backwards thru seq from i=len(seq)-1 to i=0
                    if seq[i] == 'A':
                        A_count += 1;
                    else:
                        A_count = 0
                    A_tracts[i] = A_count
                wg_index[gene] = A_tracts
                gene_len[gene] = len(seq) # for error-checking later
                gene_num += 1
            gene = mobj.group(1)
            seq=""
        else:
            seq += line
        if gene_num > 1000000:
            print "Reached gene number limit. Aborting index\n"
            break
# do last gene:
A_count = 0
A_tracts = [0] * len(seq)
for i in range(len(seq)-1, -1, -1):
    if seq[i] == 'A':
        A_count += 1;
    else:
        A_count = 0
    A_tracts[i] = A_count
wg_index[gene] = A_tracts


########################
# output of sam2table_v2 is:
#             @splut = samtable[0,2,3,5], $readlen, $cig, $ciglen, $secondary, $NH, $initial_nonmatching  (although samtable[3] used to be altered in old version)
# ie readname, transcriptname, startpos, orig_cig, orig_readlen
# want to convert to zero-based e.g. 'A' in GTAG -> (start=3; end=4; length=1). See https://www.biostars.org/p/84686/
########################

header = True # <----- HARDCODED
#header=False  # <----- HARDCODED
for fn_in in tbl_fns:
    print fn_in
    with gzip.open(fn_in,'r') as tbl:
        with gzip.open(tbl_fns[fn_in], 'w') as fout:
            if header:
                print >> fout, "# Columns: 1) gene name, 2) start, 3) end, 4) readname, 5) # hits for read,  6) strand, 7 # 5'nts not aligned, 8) # 'A's 3' of FP 9) Primary alignment? (T/F)"
            count = 0
            for aln in tbl:
                line = aln.rstrip('\n')
                pA = "NA"
                splut = line.split('\t')
                end = int(splut[2]) + int(splut[6]) 
                gene = splut[1]
                if splut[7] == 'pri':
                    splut[7] = 'T'
                else:
                    splut[7] = 'F'
                try:
                    pA = str(wg_index[gene][end-1])
                except:
                    print >> sys.stderr, "Warning: failed to find pA tract info for gene "+ gene + " position: " + str(end)
                    print >> sys.stderr, "read:" , line
                    if gene in wg_index:
                        print >> sys.stderr, "... although this gene name is in the reference."
                        print >> sys.stderr, "read start: " + str(splut[2])
                        print >> sys.stderr, "read length:" + str(splut[4])
                        print >> sys.stderr, "Next base after read end on flanked seq (zero-based):" + str(end)
                        print >> sys.stderr, "gene sequence length:" + str(gene_len[gene])
                        if gene_len[gene] == end:
                            print >> sys.stderr, "This appears to be because read mapped to very end of template sequence."
                        if gene_len[gene] < end:
                            print >> sys.stderr, "This appears to be because read alignment overhangs the end of template sequence."
                        try:
                            pA = str(wg_index[gene][gene_len[gene]-1])  # "gene_len[gene]-1" should be the index of the last element in the array in wg_index[gene]
                            print >> sys.stderr, "End pos minus one has pA.........................................", pA
                        except:
                            print >>sys.stderr, "Gene has not information on tract near end !!!!!!!!! "
                splut.extend([pA, '+', end])
                # bed format: 1-refseq_name 2-start 3-end 4-featname 5-score (I will add Num Hits / NH)  6-strand
                # custom fields: 7-leading Ss 8-pA tract len 9-primary/secondary alignment << later?? 10-pos-class 11-T/F(within_outer_limits?) 12-T/F(within_Nagalakshmi_limits) >>
                # Feature coords will be zero-based and exclusive of the 'end' coord; ATG is 1000-1003 (0-based; exclusive coords, i.e. the 'A' is the 1001st nucleotide)
                bedline=[ splut[i] for i in [1,2,12,0,8,11,9,10,7] ]
                if int(bedline[2]) == 1001:
                    if not pA=='1':
                        print >> sys.stderr, "Warning: 1nt upstream of ATG didn't give pA of 1!!! line:\n"
                        print '---'+str(len(splut))+'---'
                        print splut
                        print 'New bed file line:'
                        print bedline
                        print cdsc
                count += 1
                print  >> fout, "\t".join(str(x) for x in bedline)
