#!/bin/sh -x
# -*- mode: sh; -*-

# mangle all files into one big happy file; have to convert chromosome to fpc
# coordinates first:
convert=fgenesh-chromo2fpc.pl # in ensembl/scripts/igi

for f in "$@" ; do
    nawk -F\t '$3=="exon" || $3 ~ "_codon" ' $f |  $convert
done
# sort  is not needed, done later anyway


