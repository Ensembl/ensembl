#!/bin/sh -x
# -*- mode: sh; -*-
# $Id$

## script to do all the merges. Change to needs, this is not something very
## general or so. Although you probably could make that using some cheeky
## symlinks

## For each source (ensembl, affymetrix, fgeneh etc.), it needs one big
## collated file (called all.gtf). These files (must) have been created first
## by all-collate.sh.

ens=ensembl/all.gtf
affy=affymetrix/all.gtf
fgenesh=fgenesh/chr_gff/all.gtf 
outdir=merges
[ -d $outdir ]  && echo "Found dir $outdir, not merging" >&2  && exit 1
mkdir $outdir
prefix=igi3

gtf_merge.pl -p $prefix $ens $affy > $outdir/ens_affy.merge  2> ens_affy.log
gtf_merge.pl -p $prefix $ens $fgenesh > $outdir/ens_fgenesh.merge 2> ens_fgenesh.log
gtf_merge.pl -p $prefix $affy $fgenesh > $outdir/affy_fgenesh.merge 2> affy_fgenesh.log
gtf_merge.pl -p $prefix $ens $affy $fgenesh > $outdir/ens_affy_fgenesh.merge 2>ens_affy_fgenesh.log

cd $outdir
for i in *.merge; do
    gzip < $i > $i.gz
done

### after this, run the statistics on these files using all-stats.sh
