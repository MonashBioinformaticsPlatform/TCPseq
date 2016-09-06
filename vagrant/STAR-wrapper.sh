#!/bin/bash
IFS=' '
REGEX1="(.*)--readFilesCommand zcat(.*)"
REGEX2="(.+--readFilesIn )([^[:space:]]+)(.+)"
inCmdStr=$*
if [[ $inCmdStr =~ $REGEX1 ]];
  then inCmdStr=${BASH_REMATCH[1]}${BASH_REMATCH[2]}
  if [[ $inCmdStr =~ $REGEX2 ]];
    then
    inCmdStr=${BASH_REMATCH[1]}' <(zcat '${BASH_REMATCH[2]}') '${BASH_REMATCH[3]};
  fi;
fi;
exec /bin/bash -c "STAR.orig $inCmdStr"


