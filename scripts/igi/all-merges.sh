#!/bin/sh -x
# -*- mode: sh; -*-
# $Id$
## script to do all the merges. Change to needs, this is not something very
## general or so.
##

ens=ensembl/all.gtf
affy=2001-01-18.gs4.oo18.affymetrix/all.gtf
fgenesh=fgenesh/chr_gff/all.gtf 

# gtf_merge.pl $ens $affy > ens_affy.merge  2> ens_affy.log
# gtf_merge.pl $ens $fgenesh > ens_fgenesh.merge 2> ens_fgenesh.log
gtf_merge.pl $affy $fgenesh > affy_fgenesh.merge 2> affy_fgenesh.log
gtf_merge.pl $ens $affy $fgenesh > ens_affy_fgenesh.merge 2>ens_affy_fgenesh.log

