#!/bin/sh -x
# -*- mode: sh; -*-
# $Id$

# script to prepare all the raw data files from all the different sources
# (ensembl, affymetrix, fgeneh etc.) that will go in to the merge.  It will
# typically call a source-specific collate.sh script in the source-specific
# subdirectories

igihome=~lijnzaad/proj/igi # change to needs

cd $igihome/ensembl
ensembl.collate.sh < ensembl_oo23_october.gtf > all.gtf 2> collate.log

cd $igihome/affymetrix
affymetrix.collate.sh  > all.gtf 2> collate.log 

cd $igihome/fgenesh
fgenesh.collate.sh chr_gff/*.sgp.gff > all.gtf 2>collate.log

#  fgenesh chromo+chromo_coords may be funny, run statistics on this:
mkdir missing-stats
cd missing-stats
fgenesh.missing-stats.sh ../collate.log  # puts stuff into missing-stats/*

# add more sources here if needed

# done
