#!/bin/env bash
# -*- mode: sh; -*-
# $Id$
# Little convenience script to dump all or nearly ensembl databases by
# chromosome, using dump_by_chr.pl. Current script is meant to be edited,
# i.e., change to needs.
#

set -x # be verbose

# . ~/.ensembl-setup # set PERL5LIB etc.; write yourself
export PERL5LIB
PERL5LIB=$HOME/head/ensembl/modules:${PERL5LIB:-}
# only Bio::EnsEMBL::DBLoader; is actually used.


chrs="chrY chr19 chr20 chr21 chr22"
# chrs="chr21"
user=ensadmin
pass=ensembl
host=ensrv3
# check=
check="-check" # practically essential when running on new release for first

# 
litedb=homo_sapiens_lite_120
core=homo_sapiens_core_120
fam=homo_sapiens_family_120
disease=homo_sapiens_disease_120
maps=homo_sapiens_maps_120
expression=homo_sapiens_expression_120
snp=homo_sapiens_snp_120
embl=homo_sapiens_embl_120
est=homo_sapiens_est_120
mouse=homo_sapiens_mouse_120

# alternatively, to do everything in one fell swoop with the -all option
# (documented in dump_by_chr.pl)

all='homo_sapiens_%s_120'

here=`pwd`


timestamp() {
    echo -n "$1 " "DATE: "; date
}

# databases=" -dumplite \
#         -core $core \
#        -family $fam \
#        -disease $disease \
#        -maps $maps \
#        -expression $expression \
#        -snp $snp \
#        -embl $embl \
#        -est $est \
#        -mouse $mouse"

databases=" -family $fam \
       -disease $disease \
       -maps $maps \
       -expression $expression \
       -snp $snp \
       -embl $embl \
       -est $est"

timestamp
for chr in $chrs; do
    timestamp start
    cd $here
    mkdir $chr || true
    cd $chr
    touch ../$chr.log # in case it doesn't exist


    echo "### starting dump " >> ../$chr.log
    perl ~/head/ensembl/scripts/dump_by_chr.pl $check \
        -litedb $litedb \
        $databases \
        -u $user -p $pass  -h $host -chr $chr  2>&1 | tee -a ../$chr.log

    timestamp end
done


# ------------------------------------------------------------------------
# the alternative: 
# dump_by_chr.pl -litedb homo_sapiens_lite_120 -dumplite  \
# -all 'homo_sapiens_%s_120' -chr chr22 -h ensrv3 -u ensadmin -p ensembl
# 
