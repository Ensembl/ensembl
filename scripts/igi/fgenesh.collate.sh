#!/bin/sh -x
# -*- mode: sh; -*-

# mangle all files into one big happy file; have to convert chromosome to fpc
# coordinates first:
convert=fgenesh-chromo2fpc.pl # in ensembl/scripts/igi

infiles=*.sgp.gff

# $sort="sort -k1,1 -k7,7 -k4,4n "

for f in "$@" ; do
    nawk -F\t '$3=="exon" || $3 ~ "_codon" ' $f |  $convert
done

