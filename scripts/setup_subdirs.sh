#!/bin/bash
TCPDIR=`pwd`
mkdir sacCer3
mkdir sacCer3/orig
mkdir sacCer3/rRNA_locus
mkdir sacCer3/rRNA_locus/STAR
mkdir sacCer3/tRNA
mkdir sacCer3/tRNA/BT2
mkdir sacCer3/tRNA/STAR
mkdir sacCer3/genbank
mkdir sacCer3/structRNAs
mkdir sacCer3/structRNAs/STAR
mkdir sacCer3/mRNAs_1000flank
mkdir sacCer3/mRNAs_1000flank/STAR
mkdir sacCer3/annotation
mkdir trimmed
mkdir mapped
mkdir tabulated_data
echo 'Set up TCP-seq analysis directory structure.'
echo Before running any command, change working directory to the TCP-seq directory by typing:
echo cd $TCPDIR
echo and load user-defined variables by typing:
echo source TCPseq_config.sh
