#!/bin/sh

if [ -z $1 ] || [ -z $2 ] || [ ! -z $3 ]; then
	echo "Plots statistics of given pattern with capture";
	echo "Usage: patternstats.sh <log-file-name> <regex-pattern>";
	echo "Example: ./pattern <log-file-name> 'history_len=(\d+)'";
	exit;
fi

scripts=`dirname $0`

echo ---Stats for regex=$2---
cat $1 | perl -e "while(<>) { if(/$2/) { print \"\$1\n\"; } }" | $scripts/txthist.pl
