#!/bin/sh -x
# $Id$

# run the vital statistics on the results of the merges. This will 
# use the *.merge files that result from running all-merges.sh

igihome=$HOME/proj/igi                  # change to needs
resultdir=$igihome/out
indir=$resultdir/merged
outdir=$resultdir/stats
mappingoutdir=$resultdir/mapping

[ ! -d $indir ] && echo "$indir: not found ">&2 && exit 1
[ -d $outdir ]  && echo "Found dir $outdir, not merging" >&2  && exit 1
[ -d $mappingoutdir ] && echo "Found $mappingoutdir, not doing stats">&2 && exit 1
mkdir $outdir
mkdir $mappingoutdir
cd $indir

for m in *.merge; do
  stats-from-merge-files.pl < $m -stats -chaining 5  \
         -igi2native $mappingoutdir/$m-i2n           \
         -native2igi $mappingoutdir/$m-n2i           \
         | remap-sources.sed > $outdir/$m.stats 2> $outdir/$m.log
done


