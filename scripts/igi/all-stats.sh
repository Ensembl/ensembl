#!/bin/sh -x
# $Id$

# run the vital statistics on the results of the merges. This will 
# use the *.merge files that result from running all-merges.sh
cd merges
outdir=../stats
mappingoutdir=../mapping

for m in *.merge; do
  stats-from-merge-files.pl < $m -stats -chaining 5  \
         -igi2native $mappingoutdir/$m-i2            \
         -native2igi $mappingoutdir/$m-n2            \
         | remap-sources.sed > $outdir/$m.stats 2> $outdir/$m.log
done


