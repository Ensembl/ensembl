#!/bin/sh -x
# $Id$

# run the vital statistics on the results of the merges. This will 
# use the *.merge files that result from running all-merges.sh

igihome=$HOME/proj/igi                  # change to needs
resultdir=$igihome/out
indir=$resultdir/merged

outdir=$resultdir/stats
mappingoutdir=$resultdir/mapping
summaryoutdir=$resultdir/summary
finaloutdir=$resultdir/final

allmerges=*.merge
fullmergeonly=ens_affy_fgenesh.merge

[ ! -d $indir ] && echo "$indir: not found ">&2 && exit 1

for d in $outdir $mappingoutdir $summaryoutdir $finaloutdir; do 
    [ -d $d ]  && echo "Found dir $d, not merging" >&2  && exit 1
    mkdir $d
done

cd $indir


for m in $allmerges ; do
    name=`basename $m .merge`
    stats-from-merge-files.pl < $m -stats -chaining 5  \
         -igi2native $mappingoutdir/$name-i2n           \
         -native2igi $mappingoutdir/$name-n2i            \
         -cluster_n 2                                     \
         -gtfsummary $summaryoutdir/$name.summary          \
       > $outdir/$name.stats 2> $outdir/$name.log
done

# now produce the 'final' gtf files, i.e. the ones that are predicted by two
# or more sources:
cd $indir
for m in $fullmergeonly; do 
    name=`basename $m .merge`
    final=$name.gtf
    gtfsummary2gtf.pl $summaryoutdir/$name.summary  < $m > $finaloutdir/$final
    gzip < $finaloutdir/$final > $finaloutdir/$final.gz
done

