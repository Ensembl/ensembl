#!/bin/sh -x
# -*- mode: sh; -*-

# concatenate everything into one big happy sorted file
chroms="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y UL NA"
echo "make sure that these are the right chromosomes:" | tee /dev/tty >&2
echo $chroms >&2
echo "" > /dev/tty

for c in $chroms; do
# untested, nawk bit used to be done at end
    find $c -name '*.affymetrix.gtf' -exec nawk -F\t '$3=="exon" || $3 ~ "_codon" ' {} \; > $c/all.collated
done

# for c in $chroms; do
#     sort -k1,1 -k7,7 -k4,4n  $c/all > $c/all.sorted
# done

sort -m -k1,1 -k7,7 -k4,4n `find . -name 'all.collated' -print` 

