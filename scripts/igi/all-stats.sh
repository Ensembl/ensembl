#!/bin/sh -x
# $Id$

# run the vital statistics on the results of the merges. This will 
# use the *.merge files that result from running all-merges.sh
outdir=stats

foreach m in *.merge; do
  ./stats-from-merge-files.pl $m | remap-sources.sed   \
         > $outdir/$m.stats 2> $outdir/$m.log
done


