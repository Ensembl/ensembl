#!/bin/sh -x
# -*- mode: sh; -*-
# $Header$

# Silly script to dump the ENS[P,G] to swissprot/trembl accno mappings from an
# ensembl database. Used by gene-descriptions.pl

mysql=mysql

outname=ens-sptr-mapping
outfile=`pwd`/$outname.dat
logfile=$outname.log

Usage="Usage: $0 -h host -u user database; produces files $outname.{dat,log}";

if [ $# -le 1 ]; then
    echo $#
    echo $Usage >&2
    exit 2
fi

(cat <<EOF
SELECT tsc.translation, tsc.gene, xdb.db_name, x.dbprimary_id
FROM transcript tsc, 
     objectXref ox,
     Xref x,
     externalDB xdb
WHERE tsc.translation = ox.ensembl_id 
  AND ox.xrefId = x.xrefId
  AND x.externalDBId = xdb.externalDBId
  AND xdb.db_name in ('SWISS-PROT', 'SPTREMBL')
order by tsc.gene asc, xdb.db_name desc, x.dbprimary_id asc
into outfile '$outfile'
EOF
) | $mysql "$@" > $logfile 2>&1


