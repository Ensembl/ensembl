#!/bin/sh -x
# -*- mode: sh; -*-

# litte script to compile statistics on missing features:

outdir=./
logfile=all.log
missing=$outdir/missing
init=$outdir/assumed-on-first
stats=$outdir/missing-stats
tmp=$outdir/x

Usage="Usage: $0 reads from $logfile, writes to files '$missing, $init, $stats"

if [ $# -ne 0  ] ; then
    echo $Usage >&2
    exit 1;
fi

awk -F';' '/initial/{print $1}' $logfile | sort -u > $init
awk -F';' '/ignor/{print $1}' $logfile | sort -u > $missing

chroms="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"


echo "# features assumed to lie on first (missing) bit of fpc_contig" > $tmp
for c in $chroms; do
    echo -n "chr$c:" >> $tmp
    awk -F';' '/initial/{print $1}' $logfile | grep chr$c | sort -u | wc -l  >> $tmp
done
echo -n "# total:" >> $tmp
    awk -F';' '/initial/{print $1}' $logfile | sort -u | wc -l  >> $tmp

echo "# missing features" >> $tmp
chroms="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"
for c in $chroms; do
    echo -n "chr$c:" >> $tmp
    awk -F';' '/ignor/{print $1}' $logfile | grep chr$c | sort -u | wc -l  >> $tmp
done
echo -n "# total" >> $tmp
awk -F';' '/ignor/{print $1}' $logfile | sort -u | wc -l >> $tmp

awk -F: '/#/ || $2 > 0' < $tmp > $stats
exit 0    
