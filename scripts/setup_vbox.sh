#!/bin/bash

# Warning: not to be executed as a script.
# The below commands are intended to be executed in the interactive shell.
# ALL FOLLOWING STEPS ASSUME CWD IS cd path/to/desired/TCPseq/directory/TCPseq/ TO START WITH.

###########################################
# STEP 1.
# set up directory structure

./scripts/setup_subdirs.sh
# Notes:
# These directories are already set up on the TCPseq VM.
# However, the directories containing datafiles that are desired to be replaced should be emptied before executing the appropriate commands below.


###########################################
# STEP 2.
# user-modifiable values, will be in 'ref_config.sh' in TCPseq dir. Source it and remember to source it again if shell session is re-started:

source ref_config.sh


###########################################
# STEP 3.
# Build rRNA index (STAR)
# get chr XII 450000-491000 (note: these coordinates only apply to SacCer3 version of genome - edit in ref_config.sh as appropriate)

source ref_config.sh
./scripts/index_rRNA_locus.sh


###########################################
# STEP 4.
# Build tRNA index (bowtie2)
# a) tRNA data downloaded from http://gtrnadb.ucsc.edu
# b) Add untemplated -CCA added (minus the 'A') using python script

source ref_config.sh
cd sacCer3/tRNA
curl $TRNA_URL > $TRNA_FILE.tar.gz
gunzip $TRNA_FILE.tar.gz
tar -xf $TRNA_FILE.tar
cat $TRNA_FILE.fa | python ../../scripts/add_CC_each_seq.py > $TRNA_FILE\_CCnr.fa
# option A: bowtie2
bowtie2-build sacCer3-tRNAs_CCnr.fa BT2/sacCer3-tRNAs_CCnr
cd ../../


###########################################
# STEP 5.
# Build miscRNA reference (STAR)

source ref_config.sh
cd sacCer3
for chr in ${GBK_CHRS[@]};
do curl $GBK_PREFIX$chr$GBK_SUFFIX.gz | gunzip > genbank/$GBK_FILE_PREFIX$chr$GBK_SUFFIX;
done

perl ../scripts/make_flanked_spliced_RNAs.pl -5p 50 -3p 50 -t snoRNA misc_RNA -dir genbank --whitelist structRNAs/transposon_RNAs.txt  -out structRNAs/structRNAs.fa
STAR $STAR_BUILD --runMode genomeGenerate --genomeDir structRNAs/STAR --genomeFastaFiles structRNAs/structRNAs.fa
mv Log.out structRNAs/STAR/
cd ../


###########################################
# STEP 6.
# Build ORF + flanking reference (STAR) excluding dubious mRNAs

cd sacCer3
perl ../scripts/make_flanked_spliced_RNAs.pl -5p 1000 -3p 1000 -t mRNA -exdub yes -dir genbank -out mRNAs_1000flank/mRNAs_1000flank.fa
STAR $STAR_BUILD --runMode genomeGenerate --genomeDir mRNAs_1000flank/STAR --genomeFastaFiles mRNAs_1000flank/mRNAs_1000flank.fa
mv Log.out mRNAs_1000flank/STAR/
cd ../


###########################################
# STEP 7.
# Get gene-level information
# a) Calculate 3' UTR and 5' UTRsfrom Nagalakshmi 2008 data
# b) Flag 5' UTR introns from Zhang Z. et al 2007
# c) Look at distance to next CDS, upstream and downstream, as an upper limit to the UTR length
# d) Manually flag some repetitive genes 'YLR044C','YHR146W', 'YOR203W' (or genes which overlap with / map to structural RNAs)

source ref_config.sh
cd sacCer3/annotation
curl $UTR5P_URL > $UTR5P_FILE
curl $UTR3P_URL> $UTR3P_FILE

# This is failing for an unknown reason to download using intermine.  Stored in repo for now
#python ../../scripts/dload_5pUTR_introns.py > FPUTRIntron.txt
cp ../../scripts/FPUTRIntron.txt .

python ../../scripts/transcript_extrem.py -f $UTR5P_FILE -t $UTR3P_FILE -g ../genbank -i FPUTRIntron.txt -b "YLR044C YHR146W YOR203W" > genetable_dist_to_nextgene.txt
cd ../../

###########################################
# TROUBLESHOOTING:
# 1. STAR genomeGenerate for small genomes:
# If you get a segmentation fault, the genome may be too small for the default setting of the --genomeSAindexNbases argument.
# "Generally, --genomeSAindexNbases  needs to be scaled with the genome length, as ~min(14,log2(ReferenceLength)/2 - 1)" - A. Dobin

