#!/bin/sh -x
# -*- mode: sh; -*-
# $Id$

# script to prepare all the raw data files from all the different sources
# (ensembl, affymetrix, fgeneh etc.) that will go in to the merge.  It will
# typically call a source-specific collate.sh script in the source-specific
# subdirectories.

# If the only thing that changes across releases are file names or directory
# names, then having the right symlinks in place will make it general enough
# See comments on symlinks below each of the data sources

igihome=$HOME/proj/igi                  # change to needs
resultdir=$igihome/out
outdir=$resultdir/collated

[ -d $resultdir ]  && echo "found $resultdir, not overwriting it" >&2 && exit 1
mkdir $resultdir
mkdir $outdir

# all=all.gtf                             # where local collated data goes to
# log=collate.log                         # where local log goes to
rename=remap-sources.sed                # to translate the source field

### EnsEMBL:
source=ensembl
indir=$igihome/$source

[ ! -d $indir  ] && echo "$indir: not found" >&2 && exit 1
ensembl.collate.sh < $indir/ensembl.gtf 2> $outdir/$source.log | $rename > $outdir/$source.all
### end Ensembl


### Affymetrix
## It has everything in chromosome directories that don't seem to
## change across releases, but warn about this:
# chroms="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"
chroms="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y UL NA"
echo "make sure that these are the right chromosomes:" | tee /dev/tty >&2
echo $chroms >&2
echo "" > /dev/tty
source=affymetrix
indir=$igihome/$source
[ ! -d $indir ] && echo "$indir: not found" >&2 && exit 1
cd $indir
affymetrix.collate.sh  $chroms 2> $outdir/$source.log   | $rename > $outdir/$source.all 
### end Affymetrix

### Softberry/Fgenesh
source=fgenesh
indir=$igihome/$source                  # can be symlink
# outdir=$resultdir/$source
[ !  -d $indir ] && echo "$indir: not found" >&2 && exit 1
cd $indir

# fgenesh.collate.sh chr_gff/*.sgp.gff 2>$source.log | $rename > $source.all
fgenesh.collate.sh chr_gff/*.sgp.gff chr_r/*.sgp.gff 2>$outdir/$source.log | $rename > $outdir/$source.all

#  fgenesh chromo+chromo_coords -> fpcctg+fpc_coords mapping may be funny, run
#  statistics on this:
missingdir=$outdir/fgenesh.missing
mkdir $missingdir
cd $missingdir
fgenesh.missing-stats.sh $outdir/$source.log  # puts stuff into missing-stats/*

### end Softberry/Fgenesh


# add more sources here if needed

# done
