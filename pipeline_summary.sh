# Summary of TCP-seq pipeline commands. 
# Warning: not to be executed as a script.
# The below commands are intended to be executed in the interactive shell.

####################
#   PREPARATION :
####################

# Set up the TCPseq VM and log in
mkdir TCP_VM
cd TCP_VM
vagrant init http://bioinformatics.erc.monash.edu/TCPseq/vm/TCPseq.box
vagrant up
mkdir input_data # make a folder for the input fastq files
vagrant ssh

##################################################################################
# STEP 144
# (optional:) Update pipeline from GitHub
# From the TCPseq VM, close the pipeline_summary.sh file (and any other open file from the TCPseq repository) then type:
git pull
# Read the corresponding readme.md file.

##################################################################################
# STEP 145
# (optional:) Download test dataset:
curl http://bioinformatics.erc.monash.edu/TCPseq/vm/input_data.tar | tar -C ../input_data -xf  -
# Note: alternatively, you must transfer your fastq files to the input_data directory you created above

###########################
#  DATA ANALYSIS PIPELINE :
###########################

# Note: all steps require returning to the /vagrant/TCPseq directory

##################################################################################
# STEP 146
# Initiate certain variables.
source TCPseq_config.sh
# Notes:
# 1) Re-source this at the start of any new shell session for TCPseq pipeline.
# 2) This will report the number of input files detected the first time it runs.
# If this number is '0', delete the file 'input_filenames_auto.txt' and
# check contents of your 'input_data' directory before trying again.
# Also check that the INPUT_FQ_SUFFIX entry in the TCPseq_config.sh script
# matches the extensions of your fastq files.


##################################################################################
# STEP 147
# a) Trim poor-quality bases (where av qual < 24 over a sliding window of 7 nt), keep only reads > 26 nt long
# b) Select only those with a p(A) tract of at least 12 consecutive 'A's. Remove the p(A) tract and everything 3' of the tract.
cd trimmed
for i in $(seq $FNUM); do
  bn=${FQ_HNDL[$i]}.fq.gz
  trimmomatic SE -threads $NCORES ../$INPUT_DATADIR/${FQ_FNS[$i]} tr_$bn SLIDINGWINDOW:7:24 MINLEN:26
  cutadapt --discard-untrimmed -a AAAAAAAAAAAA -e 0.0 -O 12 tr_$bn | gzip > pA_tr_$bn
done


##################################################################################
# STEP 148.
# a) Salvage reads (~1%) having 5' adapter sequences by trimming them away.
# b) Keep only reads whose remaining sequence is 14 nt long or more
cd trimmed
for i in $(seq $FNUM); do
  bn=${FQ_HNDL[$i]}.fq.gz
  cutadapt -g TCAGAGTTCTACAGTCCGACGATC -e 0.0 -O 11 pA_tr_$bn | gzip > fpc_pA_tr_$bn
  trimmomatic SE -threads $NCORES fpc_pA_tr_$bn m14_fpc_pA_tr_$bn MINLEN:14
done


##################################################################################
# STEP 149.
# Filter out reads that map to -
# a) rRNA
# b) tRNA
# c) other structural RNAs

./scripts/stepwise_subtractive_mapping.sh


##################################################################################
# STEP 150.
# Map to ORFs + 1kb flanking

cd mapped
for i in $(seq $FNUM); do
  bn=${FQ_HNDL[$i]}
  STAR --genomeDir $MRI_S --runThreadN $NCORES --outFilterMultimapNmax 20000  --readFilesIn <(zcat noRTS_$bn.fq.gz) \
    --outFileNamePrefix noRTS_v_orf1000exdub_$bn  2> noRTS_v_orf1000exdub_$bn.txt
done


#################################################################################################################
# STEP 151.
# Convert SAM format to a table using custom Perl script sam2table_v2.pl.

cd tabulated_data
for i in $(seq $FNUM); do
  bn=${FQ_HNDL[$i]}
  cat ../mapped/noRTS_v_orf1000exdub_$bn\Aligned.out.sam | perl ../scripts/sam2table_v2.pl -min 17 \
    -out tbl1_check_$bn.sam | gzip > tbl1_$bn\_min17.txt.gz ;
done
# Notes:
# This will infer FP length from the CIGAR field according to rules in SI Fig 1.
# outputs a table w/ columns:  1-read 2-gene 3-start 4-orig_cig 5-seq_len 6-cig 7-cig_len 8-prisec 9-NH 10-leadingSs


###################################################################################################################
# STEP 152.
# i) removes alignments < 16 nt long
# ii) Add template-specified p(A) tract information for reads ending at a p(A) tract
# iii) Rearrange columns to be consistent with zero-based bed format

cd tabulated_data
python ../scripts/index_A_tracts.py -p tbl1_ -h $FQ_HNDL_STRNG -s _min17.txt.gz -m '../sacCer3/mRNAs_1000flank/mRNAs_1000flank.fa' -f 1000


###################################################################################################################
# STEP 153.
# Tweak a few things in the data prior to visualization
# a) For ambiguous 3' ends append a randomly selected 3' coordinate from all possible 3' ends
# b) Classify FPs according to location on transcript (5' UTR / start codon / ORF / Stop codon / 3' UTR)
# c) 'rgn_dens' column = Calculate background FP density of each of these regions, per-transcript, per-sample,
#                       EXCLUDING reads with the same start and end coords as the read in this row. Region size
#                       of 'start codon' and 'end codon' are 6 nt. 
# d) 'st_en_ratio' column = joint start/end frequency (per-transcript, per-sample) over 'rgn_dens'
# e) merges with gene information, UTR information
# f) Flags FPs within 30nt of downstream or upstream ORFs

cd tabulated_data
for i in $(seq $FNUM); do
  bn=${FQ_HNDL[$i]}
  Rscript ../scripts/Process_bed_for_plotting.R tbl2_$bn.bed.gz ../sacCer3/annotation/genetable_dist_to_nextgene.txt tbl3_$bn.bed.gz 1000
done
# Notes:
# 1) Outputs bed table as above, + new columns.
# 2) Zero-based bed format with 1kb flaning will place the start codon at position 1001.
#    e.g. the 'ATG' trinucleotide would be reported as 1001..1004. 


###################################################################################################################
# STEP 155.
# (Optional) Get gene-level metrics

cd tabulated_data
 Rscript ../scripts/gene_level_metrics.R ../../input_data/input_filenames_manual.txt glm_out_
# Notes:
# To filter by logical columns, add 3rd argument containing column names (preceding NOT inverts the flags) separated by commas (no spaces), e.g.:
# Rscript ../scripts/gene_level_metrics.R ../../input_data/input_filenames_manual.txt glm_out_ NOTmanually_flagged,within_Nagalakshmi_transcript
#
# UPDATE 30/1/17:
# For large files the above R process might be killed or error out. An equivalent Python script which will work with large files is now available (requires PANDAS). e.g.:
# python ../scripts/gene_level_metrics.py ../../input_data/input_filenames_manual.txt glm_out_ NOTmanually_flagged,within_Nagalakshmi_transcript

###################################################################################################################
# STEP 156. (optional)
# (Optional) Do gene-wise bootstrapping (100 resamples)
# Requires python 2.7 with numpy; pandas; docopt. Generates 3 files (start/length joint frequency; bootstrapped start marginal density; bootstrapped end marginal density) 

mkdir mtgn
for i in $(seq $FNUM); do
  bn=${FQ_HNDL[$i]}
  python scripts/metagene_agg.py -b tabulated_data/tbl3_$bn.bed.gz  -g sacCer3/annotation/genetable_dist_to_nextgene.txt -o mtgn/$bn -p 100
done

###################################################################################################################
# STEP 157. (optional)
# Make metagene plots (marginal densities etc.) from bootstrapped data generated in Step 156 above
# Requires R 3.3.2 with plyr; dplyr; ggplot2; docopt
# The call to scripts/dens_plt.R takes plotting arguments. For details, enter: Rscript scripts/dens_plt.R -h 

for i in $(seq $FNUM); do
  bn=${FQ_HNDL[$i]}
  Rscript scripts/dens_plt.R -i mtgn/$bn -l -60 -r 40
done

