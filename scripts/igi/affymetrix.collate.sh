#!/bin/sh -x
# -*- mode: sh; -*-

for chr in "$@"; do
    find $chr -name '*.affymetrix.gtf' -exec \
      nawk -F\t '$3=="exon" || $3 ~ "_codon" ' {} \;
done 
# done | sort -k1,1 -k7,7 -k4,4n : not needed, sort will take place later on
# anyway

