#!/bin/sh -x
# -*- mode: sh; -*-

# $Id$

# simple script to put all peptides files of the various sources in order.
# It collate the peptide files in case they're not one big file, and 
# makes the geneid appear as the first word of a '>...' line
# (this is only for ensembl, rest already has it like that).

# It's all done in this one script (unlike all-collates.sh, which does something
# similar for the raw gtf's, and which calls other scripts)

### sources="ensembl affymetrix fgenesh"

source=ensembl
cd $source
origfile=$source.pep.orig
pepfile=$source.pep
[ ! -f $origfile ] && echo "$origfile: Not found">&2 && exit 1
[ -f $pepfile ] && echo "found $pepfile; not overwriting">&2 && exit 1
sed '/^>/s/>\(ENSP[0-9]*\) Gene:\(ENSG[0-9]*\)\(.*$\)/>\2 \1\3/' < $origfile \
 | sed '/^>/s/^>/>ENS:/' > $pepfile    >$pepfile  || exit 2

source=affymetrix
cd $source
pepfile=$source.pep
[ -f $pepfile ] && echo "found $pepfile; not overwriting">&2 && exit 1

chroms="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y UL NA"
echo "make sure that these are the right chromosomes:" | tee /dev/tty >&2
echo $chroms >&2
echo "" > /dev/tty
for chr in "$@"; do
    find . $chr -name '*.affymetrix.aa' -exec cat {} \;
done  | sed '/^>/s/>\(.*\)\.\([0-9][0-9]*\)$/>AFFY:\1 \1.\2/' > $pepfile
# the pep file has the transcript id, which is gene_id + number. Replace it
# to make it easier for gtfsummary2pep

source=fgenesh
cd $source
pepfile=$source.pep

[ -f $pepfile ] && echo "found $pepfile; not overwriting">&2 && exit 1
for i in   chr_pro/*.sgp.pro chr_r_pro/*.sgp.pro; do
    cat $i
done  | sed '/^>C/s/^>FGENH:C/>S.C/' > $pepfile
# (replaces ">C1000004 chr1  ..."  with ">S.C1000004 chr1  ..." everywhere;
# the gtf files have 'S.C...', the pep files 'C....'. Easiest to change here.



