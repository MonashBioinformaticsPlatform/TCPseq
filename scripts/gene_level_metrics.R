#########################################
#### Prepare FILES for Shiny: ###########
#########################################

# call from Bash e.g from tabulated_data/:
# Rscript ../scripts/gene_level_metrics.R ../../input_data/input_filenames_manual.txt 
BEDFILES_PREFIX='tbl3_'
BEDFILES_SUFFIX='.bed.gz'
FILT_ARG=''

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2){
  print( "Not enough space-separated arguments given. Exiting")
  stop()
}
SAMPLES_FN = args[1]
OUT_PREFIX=args[2]
if (length(args) > 2){
  FILT_ARG = args[3]
}
if (length(args) > 3){
  BEDFILES_PREFIX = args[4]
}
if (length(args) > 4){
  BEDFILES_SUFFIX = args[5]
}

library('plyr')
library('reshape2')

filter_df = function(df, filt_cols){
  filt_all=T
  for (cn in filt_cols){
    cn2 = sub(x = cn, pattern = '^NOT', replacement = '')
    if (cn2 %in% colnames(df)){
      rowfilt = as.logical(df[,cn2])
      if (cn2 != cn ){  # invert filter as user specified a '!'
        rowfilt = !rowfilt
      }
      
    } else {
      print (paste0("Warning: column ", cn2, " not found in bed file. Ignoring."))
      rowfilt = T
    }
    filt_all = filt_all & rowfilt
  }
  return(df[filt_all,])
}

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< END FUNCTION DEFS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# for manual local execution only:
if (FALSE){
  setwd("~/VMs/bioinf/shared_data/projects/sarcher/TCPseq/tabulated_data")
  SAMPLES_FN = '../../input_data/input_filenames_manual.txt'
#  EXTREM_FN = 'sacCer3/annotation/genetable_dist_to_nextgene.txt'
  BEDFILES_PREFIX = 'tbl3_'
  BEDFILES_SUFFIX = '.bed.gz'
  OUT_PREFIX = 'genelevel_metrics_'
  FILT_ARG='!spike_flag,!manually_flagged,withinin_Nagalakshmi_transcript'
}
user_filters = strsplit(FILT_ARG, split = ',')[[1]]

# get sample info
samples2 = read.table(file = SAMPLES_FN, header = TRUE, quote = '', sep = "\t", comment.char = '')
# get gene types
#extrem = read.table(file = EXTREM_FN,
#                    quote = '', sep = "\t", header=TRUE, comment.char = '')

print ('Samples in table:')
for (samp in unique(samples2$sample)){
  print (samp)
  first=T
  for (rw in which( samples2$sample == samp )){
      tempbed = read.table(file=gzfile(paste0(BEDFILES_PREFIX, samples2$file_uniq[rw], BEDFILES_SUFFIX)), header = TRUE, quote = '', sep = "\t", comment.char = '')
      tempbed$fraction = samples2$fraction[rw]
      if (length(user_filters)>0){
        tempbed=filter_df(tempbed, user_filters)
      }
      if(first){
        all_samp_beds = tempbed 
        first=F
      } else {
        all_samp_beds = rbind(all_samp_beds, tempbed)  
      }
  }
  # do in-sample aggregating
  pc_g_f = count(all_samp_beds, c("gene", "fraction", "pos_class"))
  pc_g_f_wide = dcast(data =  pc_g_f, fun.aggregate = sum, value.var = 'freq', formula = gene ~ fraction + pos_class)
  print (paste0('Summary of ', samp, ":"))
  print(summary(pc_g_f_wide))
  write.table(pc_g_f_wide, file = paste0(OUT_PREFIX, samp, '.txt'), quote = F, sep = '\t', col.names = T, row.names = F)             
}


                            