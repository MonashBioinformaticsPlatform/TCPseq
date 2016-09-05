#!/bin/bash
source ref_config.sh

##### EDIT THE VALUES BELOW (DO NOT INTRODUCE SPACES) ##### 

export INPUT_FQ_SUFFIX=.fq.gz
export INPUT_DATADIR=../input_data
NCORES=2
SMPL_TABLE=input_filenames_auto.txt

#####  DO NOT CHANGE BELOW THIS LINE UNLESS YOU KNOW WHAT YOU ARE DOING  ######
if [ ! -r input_filenames_auto.txt ]
then
  python scripts/unique_fn_parts.py -s $INPUT_FQ_SUFFIX -d $INPUT_DATADIR -o input_filenames_auto.txt
fi
FQ_FNS=( $(cut -d $'\t' -f 1 $SMPL_TABLE) )
FQ_HNDL=( $(cut -d $'\t' -f 2 $SMPL_TABLE) )
FQ_HNDL_STRNG=${FQ_HNDL[1]}
for i in "${FQ_HNDL[@]:2}"; do
   FQ_HNDL_STRNG+=:::$i
done
FNUM=$((${#FQ_FNS[@]}-1))

# iterate thru as so:
#for i in $(seq $FNUM)
#do
#  echo ${FQ_FNS[$i]} ${FQ_HNDL[$i]} 
#done