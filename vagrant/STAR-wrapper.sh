#!/bin/bash
function join_by { local IFS="$1"; shift; echo "$*"; }
INCMDA=( "$@" )
REGEX1="(.+)--readFilesCommand zcat(.+)"
REGEX2="(.+--readFilesIn )([^[:space:]]+)(.+)"
inCmdStr=`join_by ' ' "${INCMDA[@]}"`
if [[ $inCmdStr =~ $REGEX1 ]];
  then inCmdStr=${BASH_REMATCH[1]}${BASH_REMATCH[2]}
  echo TRUE1
  echo $inCmdStr 
  if [[ $inCmdStr =~ $REGEX2 ]];
    then echo TRUE2
    inCmdStr=${BASH_REMATCH[1]}' <(zcat '${BASH_REMATCH[2]}') '${BASH_REMATCH[3]};
  fi;
fi;

/bin/bash -c "STAR.orig $inCmdStr"


