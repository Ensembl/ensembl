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

igihome=~lijnzaad/proj/igi # change to needs

cd $igihome/ensembl
ensembl.collate.sh < ensembl.gtf > all.gtf 2> collate.log
# ensembl.gft is a symlink to e.g. ensembl_oo23_october.gtf

cd $igihome/affymetrix
affymetrix.collate.sh  > all.gtf 2> collate.log 
# dir affymetrix is a symlink to e.g. 2001-02-11.gs5.oo23.affymetrix

cd $igihome/fgenesh
fgenesh.collate.sh chr_gff/*.sgp.gff > all.gtf 2>collate.log
# fgenesh is prolly als a symlink, don't know what to though :-) 

#  fgenesh chromo+chromo_coords -> fpcctg+fpc_coords mapping may be funny, run
#  statistics on this:
mkdir missing-stats
cd missing-stats
fgenesh.missing-stats.sh ../collate.log  # puts stuff into missing-stats/*


# add more sources here if needed

# done
