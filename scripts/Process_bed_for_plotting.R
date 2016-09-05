library('plyr')

# example call from BASH:
# cd tabulated_data
# Rscript ../Process_bed_for_plotting.R tbl2_0.bed.gz ../sacCer3/annotation/genetable_dist_to_nextgene.txt tbl3_0.bed.gz 1000

args <- commandArgs(trailingOnly = TRUE)
geneinfo = read.table(file = args[2], quote = '', sep = "\t", header = T)
bed = read.table(file = gzfile(args[1]), quote = '', sep = "\t",
                      colClasses = c('factor','integer','integer','character','integer','factor','integer','integer','logical'))
outfile=args[3]
if(is.na(args[4])){
  flanking = 1000
} else {
  flanking = as.integer(args[4]) 
}

# for manual testing:
if (FALSE) {
  outfile = "./bed191115_pA_gsums.txt.gz"
  out_geneinfo = "./genewise_191115_pA_gsums.txt"
  setwd('/Users/sarcher/VMs/bioinf/shared_data/projects/sarcher/TCPseq/tabulated_data')
  flanking = 1000
  bed = read.table(file = gzfile('tbl2_0.bed.gz'), quote = '', sep = "\t",
                   colClasses = c('factor','integer','integer','character','integer','factor','integer','integer','logical'))
  geneinfo = read.table(file = '../sacCer3/annotation/genetable_dist_to_nextgene.txt', quote = '', sep = "\t", header = T)
  #gene name, 2) start, 3) end, 4) readname, 5) # hits for read,  6) strand, 7 # 5'nts not aligned, 8) # 'A's 3' of FP 9) Primary alignment? (T/F)
}

colnames(bed) = c('gene','start','end','readname','NH', 'strand', 'LeadingSs','pA_tract','primary')
geneinfo = geneinfo[,c('gene','max_5p', 'max_3p', 'orf_len', 'real_5p', 'real_3p', 'intron_5pUTR', 'manually_flagged')]
# sort by transcript / pos
bed = bed[ order(bed$gene, bed$start),]
bed$min_len = bed$end - bed$start
bed$pA_tract[is.na(bed$pA_tract)] = 0
bed$max_len = bed$min_len + bed$pA_tract
set.seed(77)
bed$prob_len = sapply(bed$pA_tract,  function(x) sample(c(0:x), 1) )
bed$prob_len = bed$prob_len + bed$min_len

bed=merge(bed, geneinfo, by='gene', all.x=T)

# 6/11/15 REDEFINE 'POS_CLASS' USING CORRECTED START-END (corrected for leadingSs i.e. the CIGAR 'bug') and PENULTIMATE codon (rather than STOP codon) for 'terminating'
bed$pos_class = factor('uncertain', levels=c('uncertain', '5pUTR', 'AUG', 'ORF', 'stop', "3pUTR"))
bed$pos_class[bed$end < (flanking + 1)] = '5pUTR'
bed$pos_class[bed$end > (flanking + 4) & bed$start < (flanking -3)] = 'AUG'
bed$pos_class[bed$end < (flanking + bed$orf_len-6)+1 & bed$start > (flanking + 0)] = 'ORF'
bed$pos_class[bed$end > (flanking + bed$orf_len-6)+4 & bed$start < (flanking + bed$orf_len-6)-3 ] = 'stop'
bed$pos_class[bed$start > (flanking + bed$orf_len-6)  ] = '3pUTR'
bed$within_Nagalakshmi_transcript = T
# for FPs in the 'scanning' or '3pUTR' zone, classify as within or outside of Nagalakshmi experimentally determined UTRs 
# (or NA if no UTR data exists) 
filt = (bed$pos_class == '5pUTR')
bed$within_Nagalakshmi_transcript[filt][(bed$end[filt] + bed$pA_tract[filt]) < (flanking - bed$real_5p[filt])] = F
filt = (bed$pos_class == '3pUTR')
bed$within_Nagalakshmi_transcript[filt][(bed$start[filt]) > (flanking + bed$orf_len[filt] + bed$real_3p[filt])] = F
# flag FPs so far away that they are overlapping neighbouring CDS
bed$within_max_transcript = bed$start > (flanking - bed$max_5p) & (bed$end + bed$pA_tract) < (flanking + bed$orf_len + bed$max_3p)

# make new char column 'start_end' for start/end combo
bed['start_end'] = paste(bed$start, bed$end, sep=':')

cnt = count(bed, c("start_end", "gene", "pos_class"))
colnames(cnt)[4] = "st_en_frq_in_pc_g"
orig_col_order=c(colnames(bed),"start_end","st_en_frq_in_pc_g", "totreads_in_pc_g", "totreads_in_g")
bed <- merge(bed, cnt, by = c("start_end", "gene", "pos_class"), all.x=TRUE)

bg = count(bed, c("gene", "pos_class"))
colnames(bg)[3] = "totreads_in_pc_g"
bed <- merge(bed, bg, by = c("gene", "pos_class"), all.x=TRUE)

bg2 = count(bed, c("gene"))
colnames(bg2)[2] = "totreads_in_g"
bed <- merge(bed, bg2, by = c("gene"), all.x=TRUE)
bed=bed[,orig_col_order]

bed['rgn_dens'] = NA

logi = (bed$pos_class == '5pUTR' & bed$within_Nagalakshmi_transcript)
bed$rgn_dens[logi] = (bed$totreads_in_pc_g[logi] - bed$st_en_frq_in_pc_g[logi]) / (bed$real_5p[logi]) # background FP density EXCLUDING the start:end freq of current row
logi = (bed$pos_class == 'ORF')
bed$rgn_dens[logi] = (bed$totreads_in_pc_g[logi] - bed$st_en_frq_in_pc_g[logi]) / (bed$orf_len[logi]) # background FP density EXCLUDING the start:end freq of current row
logi = (bed$pos_class == 'AUG')
bed$rgn_dens[logi] = (bed$totreads_in_pc_g[logi] - bed$st_en_frq_in_pc_g[logi]) / (6)  # Not entirely mathematically correct.... but if FP length is 19, start coord is then within -15 : -3  = a 12 nt window
logi = (bed$pos_class == 'stop')
bed$rgn_dens[logi] = (bed$totreads_in_pc_g[logi] - bed$st_en_frq_in_pc_g[logi]) / (6)  # Not entirely mathematically correct either

bed['st_en_ratio'] = bed$st_en_frq_in_pc_g / bed$rgn_dens  # this start/end combo divided by FP density in region
bed['spike_flag'] = F   # decide these thresholds in graphing program
bed$spike_flag[ bed$st_en_ratio > 20 & bed$st_en_frq_in_pc_g_s > 20 ] = T  
write.table(x = bed, file = gzfile(outfile), sep = "\t", 
            col.names = TRUE, quote = FALSE, row.names = FALSE )

#make_genewise?
 
