#########################################
#### Prepare FILES for Shiny: ###########
#########################################

# call from Bash e.g from tabulated_data/:
# python ../scripts/gene_level_metrics.py ../../input_data/input_filenames_manual.txt 
import sys
import gzip
import numpy as np
import pandas as pd

BEDFILES_PREFIX='tbl3_'
BEDFILES_SUFFIX='.bed.gz'
FILT_ARG=''
POS_CLASSES=['3pUTR','5pUTR','AUG','ORF','stop','uncertain']
# for manual local execution only:
if False:
  SAMPLES_FN = '../../input_data/input_filenames_manual.txt'
#  EXTREM_FN = 'sacCer3/annotation/genetable_dist_to_nextgene.txt'
  BEDFILES_PREFIX = 'tbl3_'
  BEDFILES_SUFFIX = '.bed.gz'
  OUT_PREFIX = 'genelevel_metrics_'
  FILT_ARG='NOTspike_flag,NOTmanually_flagged,within_Nagalakshmi_transcript,Blah'

args = sys.argv
if len(args) < 2:
    print "Not enough space-separated arguments given. Exiting"
    exit()
SAMPLES_FN = args[1]
OUT_PREFIX=args[2]
if len(args) > 3:
    FILT_ARG = args[3]
if len(args) > 4:
    BEDFILES_PREFIX = args[4]
if len(args) > 5:
    BEDFILES_SUFFIX = args[5]

# parse user filters
user_filters = FILT_ARG.split(',')
logi = [x.startswith('NOT') for x in user_filters]
F_ind = np.where(np.array(logi))[0] # Returns indices of filters in 'user_filters' on which to filter for FALSE
user_filters_F = [user_filters[x][3:] for x in F_ind] # Returns colnames in BED file on which to filter for FALSE
T_ind = np.where(np.array([not x for x in logi]))[0] # Returns indices of filters in 'user_filters' on which to filter for TRUE
user_filters_T = [user_filters[x] for x in T_ind] # Returns colnames in BED file on which to filter for TRUE

# import samples info table
samples2 = pd.read_table(SAMPLES_FN)
# ensure all cols are interpreted as strings
samples2 = samples2.astype('str')
samples = pd.Series.unique(samples2['sample'])

for smp in samples:
    print 'Processing sample: ' + smp
    first=True
    gene_dict=dict()
    smp_rws=np.where(samples2['sample']==smp)[0]
    smp_frct=pd.Series.unique(samples2['fraction'][smp_rws])# get sub-fractions (80S, SSU etc) from a sample
    smp_frct.sort()
    colnames=list()
    for frct in smp_frct:
        colnames = colnames + [frct+'_'+pc for pc in POS_CLASSES]
    for rw in smp_rws: 
        fin_nm = BEDFILES_PREFIX + samples2['file_uniq'][rw] + BEDFILES_SUFFIX
        print 'Processing file: ' + fin_nm
        frct = samples2['fraction'][rw]
        col_offset=min([i for i, x in enumerate(colnames) if x.startswith(frct)]) # find first column starting with 'frct_'
        # we will only be adding to these columns from the current file: col_offset:col_offset+len(POS_CLASSES)
        with gzip.open(fin_nm, 'r') as fin:
            excluded=0
            included=0
            first=True
            for line in fin:
                splut = line.rstrip('\n').split('\t')
                if first:
                    #0-gene	1-start	2-end	3-readname	4-NH	5-strand	6-LeadingSs	7-pA_tract	8-primary	9-min_len	10-max_len	11-prob_len	12-max_5p	13-max_3p	14-orf_len	15-real_5p	16-real_3p	17-intron_5pUTR	18-manually_flagged	19-pos_class	
                    #20-within_Nagalakshmi_transcript	21-within_max_transcript	22-start_end	23???-start_end.1	24-st_en_frq_in_pc_g	25-totreads_in_pc_g	26-totreads_in_g	27-rgn_dens	28-st_en_ratio	29spike_flag
                    first=False
                    T_logi=[x in user_filters_T for x in splut] # =logical vector of same len as splut
                    T_filt_ind = np.where(np.array(T_logi))[0]
                    F_logi=[x in user_filters_F for x in splut] # =logical vector of same len as splut
                    F_filt_ind = np.where(np.array(F_logi))[0]
                    continue
                T_filtlist = [splut[i] for i in T_filt_ind]
                F_filtlist = [splut[i] for i in F_filt_ind]
                if any([x != 'TRUE' for x in T_filtlist]):
                    excluded +=1
                    continue
                if any([x != 'FALSE' for x in F_filtlist]):
                    excluded+=1
                    continue
                included+=1
                gene=splut[0]
                if not gene in gene_dict:
                    gene_dict[gene] = np.zeros(len(colnames), dtype=np.int)
                posclass = splut[19]
                col_offset2=min([i for i, x in enumerate(POS_CLASSES) if x == posclass])
                gene_dict[gene][col_offset+col_offset2] += 1
            print 'Excluded ' + str(excluded)+ ' of ' + str(included+excluded) + ' reads.'
    with open(OUT_PREFIX+smp+ '.txt', 'w') as fout:
        print >> fout, 'gene\t'+'\t'.join(colnames)
        allgns = [x for x in gene_dict]
        allgns.sort()
        for gene in allgns:
            print >> fout, gene+'\t'+'\t'.join([str(x) for x in gene_dict[gene].tolist()])


