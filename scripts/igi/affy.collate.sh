#!/bin/sh -x
# -*- mode: sh; -*-

# used to concatenate everything into one big happy sorted file
chroms="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"
for c in $chroms; do
    find $c -name '*.affymetrix.gtf' -exec cat {} \; > $c/all
done

for c in $chroms; do
    sort -k1,1 -k7,7 -k4,4n  $c/all > $c/all.sorted
done

all=./all.gtf

# sort it and throw everything that's not an exon or start/stop codon:
sort -m -k1,1 -k7,7 -k4,4n `find . -name 'all.sorted' -print`  | nawk -F\t '$3=="exon" || $3 ~ "_codon" '  > $all

