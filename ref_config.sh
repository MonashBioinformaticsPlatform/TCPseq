#!/bin/bash

##### EDIT THE VALUES BELOW (DO NOT INTRODUCE SPACES) ##### 

NCORES=2
FA_URL=ftp://ftp.ensembl.org/pub/release-85/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz
FA_FILE=Saccharomyces_cerevisiae.R64-1-1.dna_sm.toplevel.fa.gz #$TCPDIR/sacCer3/orig/
RRNA_URL='XII dna_sm:chromosome chromosome:R64-1-1:XII:1:1078177:1'
RRNA_START=491000
RRNA_LEN=41001
RRNA_FILE=SC3_chrXII_450000_491000.fa #$TCPDIR/sacCer3/rRNA_locus/
TRNA_URL=http://gtrnadb.ucsc.edu/genomes/eukaryota/Scere3/sacCer3-tRNAs.tar.gz
TRNA_FILE=sacCer3-tRNAs  #$TCPDIR/sacCer3/tRNA/
GBK_CHRS=(I II III IV V VI VII VIII IX X XI XII XIII XIV XV XVI Mito)
GBK_PREFIX=http://ftp.ensembl.org/pub/release-71/genbank/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.EF4.71.chromosome.
GBK_SUFFIX=.dat
GBK_FILE_PREFIX=Saccharomyces_cerevisiae.EF4.71.chromosome. #$TCPDIR/sacCer3/genbank/
UTR5P_URL=http://downloads.yeastgenome.org/published_datasets/Nagalakshmi_2008_PMID_18451266/track_files/Nagalakshmi_2008_5UTRs_V64.bed 
UTR5P_FILE=Nagal_5pUTRs.bed # $TCPDIR/sacCer3/annotation/
UTR3P_URL=http://downloads.yeastgenome.org/published_datasets/Nagalakshmi_2008_PMID_18451266/track_files/Nagalakshmi_2008_3UTRs_V64.bed
UTR3P_FILE=Nagal_3pUTRs.bed # $TCPDIR/sacCer3/annotation/
STAR_BUILD="--genomeChrBinNbits 8 --genomeSAindexNbases 10  --genomeSAsparseD 2"



############ DO NOT EDIT BELOW THIS LINE UNLESS YOU KNOW WHAT YOU ARE DOING ############
RRI_S=../sacCer3/rRNA_locus/STAR
TRI_B=../sacCer3/tRNA/BT2/sacCer3-tRNAs_CCnr
SRI_S=../sacCer3/structRNAs/STAR
MRI_S=../sacCer3/mRNAs_1000flank/STAR/



