#!/bin/sh -x
# -*- mode: sh; -*-
# $Id$

## script to do all the merges. Change to needs, this is not something very
## general or so. Although you probably could make that using some cheeky
## symlinks

## For each source (ensembl, affymetrix, fgeneh etc.), it needs one big
## collated file (called all.gtf). These files (must) have been created first
## by all-collate.sh.

igihome=$HOME/proj/igi                  # change to needs
resultdir=$igihome/out
outdir=$resultdir/merged
indir=$resultdir/collated

prefix=igi3

[ -d $outdir ]  && echo "Found dir $outdir, not merging" >&2  && exit 1
mkdir $outdir
[ ! -d $indir ] && echo "$indir: not found" >&2 && exit 1

ens=ensembl.all
affy=affymetrix.all
fgenesh=fgenesh.all

cd $outdir
gtf_merge.pl -p $prefix $indir/$ens $indir/$affy > ens_affy.merge  2> ens_affy.log
gtf_merge.pl -p $prefix $indir/$ens $indir/$fgenesh > ens_fgenesh.merge 2> ens_fgenesh.log
gtf_merge.pl -p $prefix $indir/$affy $indir/$fgenesh > affy_fgenesh.merge 2> affy_fgenesh.log
gtf_merge.pl -p $prefix $indir/$ens $indir/$affy $indir/$fgenesh > ens_affy_fgenesh.merge 2>ens_affy_fgenesh.log

cd $outdir
for i in *.merge; do
    gzip < $i > $i.gz
done

### after this, run the statistics on these files using all-stats.sh
