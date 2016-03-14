#!/bin/sh

if [ -z $1 ]; then
	echo "stability statistics plotter";
	echo "Usage: stabstats.sh <log-file-name>";
	exit;
fi

scripts=`dirname $0`

echo ---Ls---
cat $1 | perl -e "while(<>) { if(/^Ls=(\d+)/) { print \"\$1\n\"; } }" | $scripts/txthist.pl
echo ---t_bcp---
cat $1 | perl -e "while(<>) { if(/^Ls.* t_bcp=(\d+)/) { print \"\$1\n\"; } }" | $scripts/txthist.pl
echo ---Lc---
cat $1 | perl -e "while(<>) { if(/^Ls.* Lc=(\d+)/) { print \"\$1\n\"; } }" | $scripts/txthist.pl
