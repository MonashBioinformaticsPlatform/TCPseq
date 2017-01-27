HAVEPANDAS = False  # will install pandas on VM in a later commit
import getopt
import sys
import numpy as np
import os
import gzip
if HAVEPANDAS:
    import pandas as pd

#CWD = '/Users/sarcher/VMs/bioinf/shared_data/projects/sarcher/vagrant_sandbox_230117/TCP_VM/TCPseq/'
#IN_BED='./tabulated_data/tbl3_169.bed.gz'
GENE_INFO='./sacCer3/annotation/genetable_dist_to_nextgene.txt'
MAX_X = 1000
MIN_X = -1000
FLANKING=1000
MAX_LEN = 100 # zero-based: len range will be from 0 to 100
ANCHOR = 'start' # choices: 'start','stop','fp', 'tp'
MAX_BS=5
OUT_FN_PREF='./MG_test'
START_SEED=1

options, remainder = getopt.getopt(sys.argv[1:], 'b:g:l:r:f:x:a:p:o:s:c:')
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
            #print splut[1]
    
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

genes_not_found=0
with gzip.open(IN_BED, 'r') as fin:
    first=True
    for line in fin:
        if first:
            first=False
            continue
        splut = line.decode().rstrip('\n').split('\t')       
        start=int(splut[1])
        end=int(splut[2])
        flen=end-start
        try:
            gn_data=allgenes_dc[splut[0]] # gn_data[0]=bs_index ; gn_data[1] = orf_len 
        except:
            genes_not_found += 1
            continue
        curr_gweight = bsweights[gn_data[0],:]
        if ANCHOR == 'stop':  # not tested for OBO errors
            start=start-gn_data[1]
            end=end-gn_data[1]
        if start-FLANKING < MAX_X and start-FLANKING > MIN_X: #<--- -1??
            if flen+1 < MAX_LEN:
                joint_agg[flen,start]+=1
            start_marg[:,start] += curr_gweight
        if end-FLANKING < MAX_X and end-FLANKING > MIN_X: #<--- -1??
            end_marg[:,end] += curr_gweight
if genes_not_found > 0:
    print "Warning: " + str(genes_not_found) + " reads could not be assigned to genes in " + GENE_INFO
    
x_axis = np.arange(MIN_X, MAX_X).reshape(1,MAX_X-MIN_X)
joint_agg=np.append(x_axis-1, joint_agg , axis=0)   # adjust x-axis by 1 to counter out-by-one error 
start_marg=np.append(x_axis-1, start_marg , axis=0)  # adjust x-axis by 1 to counter out-by-one error
end_marg=np.append(x_axis-2, end_marg , axis=0)    # adjust x-axis by 2 to counter 2x out-by-one errors with fragment end coords

np.savetxt(fname=OUT_FN_PREF+'_jnt.txt', X=np.transpose(joint_agg),fmt='%1d', delimiter='\t')
np.savetxt(fname=OUT_FN_PREF+'_strt.txt', X=np.transpose(start_marg),fmt='%1d', delimiter='\t')
np.savetxt(fname=OUT_FN_PREF+'_end.txt', X=np.transpose(end_marg),fmt='%1d', delimiter='\t')

