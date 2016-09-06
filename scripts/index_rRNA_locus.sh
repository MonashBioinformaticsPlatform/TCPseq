#!/bin/bash

source ref_config.sh
cd sacCer3
curl $FA_URL > orig/$FA_FILE
echo '>chrXII_450000_491000' > rRNA_locus/$RRNA_FILE
zcat orig/$FA_FILE | \
    awk '/'"$RRNA_URL"' REF/{flag=1;next}/>/{flag=0}flag'  | \
    tr -d '\n' | \
    head -c $RRNA_START | \
    tail -c $RRNA_LEN | \
    fold -w60 >> \
    rRNA_locus/$RRNA_FILE
echo "" >> rRNA_locus/$RRNA_FILE
STAR $STAR_BUILD --runMode genomeGenerate --genomeDir rRNA_locus/STAR --genomeFastaFiles rRNA_locus/$RRNA_FILE
mv Log.out rRNA_locus/STAR

