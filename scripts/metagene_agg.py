HAVEPANDAS = False  # will install pandas on VM in a later commit
import getopt
import sys
import numpy as np
import os
import gzip
import time
start_time=time.time()
if HAVEPANDAS:
    import pandas as pd
    
def coerce_num(s):
    try:
        return int(s)
    except ValueError:
        return float(s)

def make_comparison(x, op, thresh, na_to_false=True):
    if na_to_false and x == 'NA':
        return(False)
    if op == '==':
        return(x == thresh)
    if op == '<':
        return(coerce_num(x) < thresh) # assumes 'thresh' has already been coerced to num!
    if op == '>':
        return(coerce_num(x) > thresh)

GENE_INFO='./sacCer3/annotation/genetable_dist_to_nextgene.txt'
MAX_X = 1000
MIN_X = -1000
FLANKING=1000
MAX_LEN = 100 # zero-based: len range will be from 0 to 100
ANCHOR = 'start' # choices: 'start','stop','fp', 'tp'
MAX_BS=5
OUT_FN_PREF='./MG_test'
START_SEED=1
uORF_FN='' #Ingolia_Table_S3_uORF.txt
MIN_SCORE=0
INSTRING='fred>10;wilma==FALSE'
INFILTS=''
options, remainder = getopt.getopt(sys.argv[1:], 'b:g:l:r:f:x:a:p:o:s:c:u:e:z:h:')
#print options

for opt, arg in options:
    if opt in ('-b'):
        IN_BED=arg
    elif opt in ('-g'):
        GENE_INFO=arg
    elif opt in ('-l'):
        MIN_X = int('-'+arg)
    elif opt in ('-r'):
        MAX_X = int(arg)  
    elif opt in ('-f'):
        FLANKING = int(arg)
    elif opt in ('-x'):
        MAX_LEN = int(arg)
    elif opt in ('-a'):
        ANCHOR = arg
    elif opt in ('-p'):
        MAX_BS = int(arg)
    elif opt in ('-o'):
        OUT_FN_PREF=arg
    elif opt in ('-s'):
        START_SEED=int(arg)
    elif opt in ('-c'):
        CWD=arg
        os.chdir(CWD)
    elif opt in ('-u'):
        uORF_FN=arg
    elif opt in ('-e'):
        MIN_SCORE = float(arg)
    elif opt in ('-z'):
        INFILTS=arg.split(',')
    elif opt in ('-h'):
        LEN_HIST_COL=arg

col_ineqs=dict()
comparison_ops=list(['==','<','>'])
for entry in INFILTS:
    for co in comparison_ops:
        if co in entry:
            splut=entry.split(co)
            if co == '==':
                col_ineqs[splut[0]] = list([co, splut[1]])
            else:
                col_ineqs[splut[0]] = list([co, coerce_num(splut[1])])
    
allgenes_dc = dict()
if HAVEPANDAS:
    genes_df = pd.read_table(GENE_INFO,sep='\t')
    ngene = len(genes_df.axes[0])
    allgenes=genes_df['gene']
    genes_df.index=genes_df['gene']
    # convert to dictionary for rapid lookup of [index, orf_len] using gene-name
    for i in range(0,ngene):
        allgenes_dc[allgenes[i]]=[i , genes_df.loc[allgenes[i], 'orf_len']]
else:
    with open(GENE_INFO, 'r') as fin:
        i=0
        for line in fin:
            splut=line.decode().rstrip('\n').split('\t')
            if splut[9]=='orf_len':
                continue
            allgenes_dc[splut[0]] = [i, int(splut[9])]
            i+=1
    ngene=i
# optional: load uORF coords.
# col0=gene name; col1=uORF start rel ORF start; col2(optional)=score (must be > MIN_SCORE)
# Note: can be multiple per gene.
gene_uORF_offsets = dict()
if uORF_FN != '':
    with open(uORF_FN, 'r') as fin:
        for line in fin:
            splut = line.rstrip('\n').split('\t')
            if len(splut) > 2 and float(splut[2]) < MIN_SCORE:
                continue
            gene=splut[0]
            if not gene in gene_uORF_offsets:
                gene_uORF_offsets[gene] = [int(splut[1])]
            else:
                gene_uORF_offsets[gene].extend([int(splut[1])])

# make gene-weights matrix with 1 column per boostrap iteration
bsweights = np.zeros(MAX_BS*ngene, dtype=np.int).reshape(ngene,MAX_BS) # rows = genes, cols = BS
for bs in range(0,MAX_BS):
    np.random.seed(seed=bs+START_SEED)
    rsamp = np.random.choice(ngene, ngene)
    bsweights[:,bs]=np.bincount(rsamp, minlength=ngene)
# add a left-hand column of just ones, which will yield the results from the non-bootstrapped data
bsweights = np.append(np.ones(ngene, dtype=np.int).reshape(ngene,1), bsweights, axis=1) 

joint_agg = np.zeros((MAX_X-MIN_X)*(MAX_LEN+1), dtype=np.int).reshape((MAX_LEN+1),(MAX_X-MIN_X)) # y-axis is flen
start_marg = np.zeros((MAX_X-MIN_X)*(MAX_BS+1), dtype=np.int).reshape(MAX_BS+1,(MAX_X-MIN_X)) # y-axis is bootstrap #
end_marg = np.zeros((MAX_X-MIN_X)*(MAX_BS+1), dtype=np.int).reshape(MAX_BS+1,(MAX_X-MIN_X)) # y-axis is bootstrap #

reads_counted=0
tot_reads=0
excluded=0
gns_19_1003 = dict()
col_ineqs_indices=dict()
len_hist_colnum=-1
print 'Start Time: ' + str(time.time()-start_time)

with gzip.open(IN_BED, 'r') as fin:
    first=True
    for line in fin:
        splut = line.decode().rstrip('\n').split('\t')
        if first:  # process column names
            first=False
            for i in range(len(splut)):
                if splut[i] in col_ineqs:
                    col_ineqs_indices[i]=col_ineqs[splut[i]]  #
                if splut[i] == LEN_HIST_COL:
                    len_hist_colnum=i
            if LEN_HIST_COL != '' and len_hist_colnum == -1:
                print "Warning: length hist column name not found in bed file. Will not produce subcategorized fragment length data."
            continue
        tot_reads+=1
        if tot_reads % 1000000 == 0:
             print  "Reads processed: " + str(tot_reads) + ". Time: "+ str(time.time()-start_time) + "  (from file "+ IN_BED + ")."
        excludethis=False
        for filtcol in col_ineqs_indices:
             if not make_comparison(x = splut[filtcol], op = col_ineqs_indices[filtcol][0], thresh = col_ineqs_indices[filtcol][1]):
                 excludethis=True
        if excludethis:
            excluded+=1
            continue
        if ANCHOR == 'uORFs':
            if not splut[0] in gene_uORF_offsets:
                continue
            gn_offsets = np.array(gene_uORF_offsets[splut[0]])
        else:
            gn_offsets = np.array([0]) # by default, already anchored to AUG + 1000 (FLANKING) 
        start=int(splut[1])
        flen=int(splut[2])-start
        try:
            gn_data=allgenes_dc[splut[0]] # gn_data[0]=bs_index ; gn_data[1] = orf_len 
        except:
            continue  # need gn_data for bootstrap weight vector of gene (not just orf_len)
        curr_gweight = bsweights[gn_data[0],:]
        if ANCHOR == 'stop':  # not yet tested for out-by-one errors
            gn_offsets = np.array([gn_data[1]])
        start=start-gn_offsets # convert start to np.array, can have >1 value for uORFs
        for st in start:
            end=flen+st
            if end-FLANKING < MAX_X and end-FLANKING > MIN_X: #<--- -1??
                end_marg[:,end] += curr_gweight
            if st-FLANKING < MAX_X and st-FLANKING > MIN_X: #<--- -1??
                start_marg[:,st] += curr_gweight
                if flen+1 < MAX_LEN:
                    joint_agg[flen,st]+=1
        reads_counted+=1

print str(reads_counted) + " reads included in output."
if reads_counted < tot_reads:
    print "Warning: " + str(tot_reads-reads_counted) + " reads (of "+ str(tot_reads) +") could not be assigned to genes using provided filters and gene information file."
    print str(excluded) + " of these were excluded by user-defined filters."
    
x_axis = np.arange(MIN_X, MAX_X).reshape(1,MAX_X-MIN_X)
joint_agg=np.append(x_axis-1, joint_agg , axis=0)   # adjust x-axis by 1 to counter out-by-one error 
start_marg=np.append(x_axis-1, start_marg , axis=0)  # adjust x-axis by 1 to counter out-by-one error
end_marg=np.append(x_axis-2, end_marg , axis=0)    # adjust x-axis by 2 to counter 2x out-by-one errors with fragment end coords

np.savetxt(fname=OUT_FN_PREF+'_jnt.txt', X=np.transpose(joint_agg),fmt='%1d', delimiter='\t')
np.savetxt(fname=OUT_FN_PREF+'_strt.txt', X=np.transpose(start_marg),fmt='%1d', delimiter='\t')
np.savetxt(fname=OUT_FN_PREF+'_end.txt', X=np.transpose(end_marg),fmt='%1d', delimiter='\t')

print "Time taken to process file " + IN_BED + ": " + str(time.time()-start_time)
