#!/bin/sh -x
# -*- mode: sh; -*-
# $Id$

# Script to merge existing gene_descriptions with those from the families
# (where known) into a new table called $new_gene_description, which is a
# complete replacement for the old table.  The script assumes that the family
# and ensembl-core database live in the same server (so you can do joins).

usage="Usage: $0 core_database family_database MUSG -h host -u user"
if [ $# -lt 4 ] ; then
    echo $usage
    exit 1
fi

# where mysql lives (or how it can be found, provided PATH is ok):
mysql=mysql; mysql_extra_flags='--batch'
# mysql=cat  # for debugging

# Database to read existing descriptions from:
core_db=$1; shift

# family database:
fam_db=$1; shift;
# check the name contains 'fam', to avoid trouble:
if echo "$fam_db" | grep -i 'fam' > /dev/null 2>&1 ; then
    : # OK
else
    echo "arg 2:  '$fam_db' does not match 'fam'" >&2
    exit 2
fi

# inside the family database, which database to read additional descriptions
# from (e.g, ENSG, ENSMUSG etc.; it's the prefix of the IDs)
prefix=$1; shift;

# check if the prefix makes sense
should_match='^ENS([A-Z]{3})?G$' # egrep pattern
if echo "$prefix" | egrep $should_match > /dev/null 2>&1 ; then
    : # OK
else
    echo "arg 2: database name to use: '$prefix' does not match $should_match" >&2
    exit 2
fi

# Name of the table having the existing descriptions (only changes when the
# schema does, really, but also useful during testing/debugging)
core_gene_desc_table='gene_description'

# the table with new descriptions to be created (in the current database):
merged_gene_desc_table='gene_description'

# now produce the SQL (with the $variables  being replaced with their values)
# and pipe this as input into mysql:
(cat <<EOF
SELECT COUNT(*) AS all_old_descriptions
FROM $core_db.$core_gene_desc_table gd
WHERE gd.description IS NOT NULL
  AND gd.description NOT IN ('', 'unknown', 'UNKNOWN');

# creating new local gene_description table with same layout+types as main one:
CREATE TABLE $merged_gene_desc_table
SELECT *
FROM  $core_db.$core_gene_desc_table gd 
WHERE gd.gene_id IS NULL;

ALTER TABLE $merged_gene_desc_table ADD PRIMARY KEY(gene_id);

# insert existing desc's (from swissprot)
INSERT INTO $merged_gene_desc_table 
  SELECT *
  FROM $core_db.$core_gene_desc_table;

# delete anything that looks remotely unknown:
DELETE FROM $merged_gene_desc_table WHERE description is null;
DELETE FROM $merged_gene_desc_table WHERE description = '';
DELETE FROM $merged_gene_desc_table WHERE description = 'unknown';
DELETE FROM $merged_gene_desc_table WHERE description = 'UNKNOWN';
DELETE FROM $merged_gene_desc_table WHERE description REGEXP '^[ \t]*$';
DELETE FROM $merged_gene_desc_table WHERE description REGEXP '^.$';

# selecting all genes from gene, this time properly as 'unknown':
INSERT INTO $merged_gene_desc_table
  SELECT id, 'unknown'
  FROM $core_db.gene;
# this will fail on the knowns, as they should; only the 'unknowns' do get
# in. This is so that all the unknowns are recognizable by the string 'unknown'.

# give stats:
SELECT COUNT(*) AS unknown_old_descriptions
FROM $merged_gene_desc_table
WHERE description = 'unknown';
  
# set up a temporary table to hold descriptions from the family database, for
# those ones that are unknown.  (This table is local to our mysql session,
# and is deleted automatically when it ends).
CREATE TEMPORARY TABLE tmp_new_descriptions
  SELECT gd.gene_id, f.description
  FROM family f, family_members fm, $merged_gene_desc_table gd
  WHERE gd.description ='unknown'
    AND fm.db_name ='$prefix'
    AND fm.db_id = gd.gene_id
    AND fm.family = f.internal_id
    AND f.description <> 'UNKNOWN';

# get rid of the unknowns:
DELETE FROM $merged_gene_desc_table where description = 'unknown';

# and insert them from the tmp (this should produce any duplicate
# warnings. But that wouldn't matter anyway). 
INSERT INTO $merged_gene_desc_table
  SELECT *
  FROM tmp_new_descriptions;

# done; give new stats:
SELECT COUNT(*) AS all_new_descriptions
FROM $merged_gene_desc_table;
EOF
) | $mysql  $mysql_extra_flags "$@" $fam_db
