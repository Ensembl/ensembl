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

echo "\n*** This script is deprecated. It has been integrated in gene-descriptions.pl ***\n";

if [ $# -le 1 ]; then
    echo $#
    echo $Usage >&2
    exit 2
fi

(cat <<EOF
SELECT tsc.translation_id, tsc.gene_id, xdb.db_name, x.dbprimary_id, ix.query_identity, ix.target_identity
FROM transcript tsc, 
     objectXref ox,
     Xref x,
     externalDB xdb,
     identityXref ix
WHERE tsc.translation_id = ox.ensembl_id 
  AND ox.xrefId = x.xrefId
  AND x.externalDBId = xdb.externalDBId
  AND xdb.db_name in ('SWISSPROT', 'SPTREMBL')
  AND ox.objectxrefId = ix.objectxrefId
order by tsc.gene_id asc, xdb.db_name desc, x.dbprimary_id asc
into outfile '$outfile'
EOF
) | $mysql "$@" > $logfile 2>&1


