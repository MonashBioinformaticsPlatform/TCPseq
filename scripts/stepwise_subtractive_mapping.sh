#!/bin/bash
source TCPseq_config.sh
cd mapped
for i in $(seq $FNUM); do
  bn=${FQ_HNDL[$i]}
  # Map to rRNA locus with STAR (both strands)
  STAR --outSAMunmapped Within --genomeDir $RRI_S --runThreadN $NCORES --outFilterMultimapNmax 20000  \
    --readFilesIn <(zcat ../trimmed/m14_fpc_pA_tr_$bn.fq.gz) --outFileNamePrefix fpc_vs_rR_$bn  2> fpc_vs_rR_$bn.txt
  # get unmapped reads
  samtools view -S -f 4 fpc_vs_rR_$bn\Aligned.out.sam | awk '{print "@"$1"\n"$10"\n+\n"$11}' | gzip > noR_$bn.fq.gz;
  # map to tRNA (stranded) using sensitive settings on Bowtie2; collect unmapped reads directly from Bowtie2 using the --un option
  bowtie2 --norc -a --very-sensitive -x $TRI_B -U noR_$bn.fq.gz --un norRtR_$bn.fq 2>$bn\_bt_err.txt | samtools view -bSh - > tR_vs_norR_$bn\.bam
  gzip norRtR_$bn.fq
  # map to structural RNAs with STAR (stranded)
  STAR --outSAMunmapped Within --genomeDir $SRI_S --runThreadN $NCORES --outFilterMultimapNmax 20000  \
    --readFilesIn <(zcat norRtR_$bn.fq.gz) --outFileNamePrefix RTS_$bn  2>miscmap_star_$bn.txt
  # get unmapped reads and also reads that mapped only in reverse orientation
  python ../scripts/rc_only_from_sam.py -i RTS_$bn\Aligned.out.sam -f ../trimmed/m14_fpc_pA_tr_$bn.fq.gz | gzip > noRTS_$bn.fq.gz
done
