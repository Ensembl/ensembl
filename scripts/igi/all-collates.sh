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

all=all.gtf                             # where local collated data goes to
log=collate.log                         # where local log goes to
rename=remap-sources.sed                # to translate the source field

### EnsEMBL:
dir=$igihome/ensembl                    # can be symlink
[ ! -d $dir ] && echo "$dir: not found" >&2 && exit 1
cd $dir
[ -f $all  ] && echo "found $all, not overwriting it ">&2 && exit 1
ensembl.collate.sh < ensembl.gtf 2> $log | $rename > $all
### end Ensembl


### Affymetrix
## It has everything in chromosome directories that don't seem to
## change across releases, but warn about this:
chroms="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"  #UL NA"
echo "make sure that these are the right chromosomes:" | tee /dev/tty >&2
echo $chroms >&2
echo "" > /dev/tty
dir=$igihome/affymetrix                    # can be symlink
[ ! -d $dir ] && echo "$dir: not found" >&2 && exit 1
cd $dir
[ -f $all  ] && echo "found $all, not overwriting it ">&2 && exit 1
affymetrix.collate.sh  $chroms 2> $log   | $rename > $all 
### end Affymetrix

### Softberry/Fgenesh
dir=$igihome/fgenesh                    # can be symlink
[ !  -d $dir ] && echo "$dir: not found" >&2 && exit 1
cd $dir
[ -f $all  ] && echo "found $all, not overwriting it ">&2 && exit 1

c2f_log=chromo2fpc.log
# fgenesh.collate.sh chr_gff/*.sgp.gff chr_r/*.sgp.gff  > $all 2>$log
fgenesh.collate.sh chr_gff/*.sgp.gff 2>$log | $rename > $all

#  fgenesh chromo+chromo_coords -> fpcctg+fpc_coords mapping may be funny, run
#  statistics on this:
mkdir missing-stats
cd missing-stats
fgenesh.missing-stats.sh ../collate.log  # puts stuff into missing-stats/*


### end Softberry/Fgenesh


# add more sources here if needed

# done
