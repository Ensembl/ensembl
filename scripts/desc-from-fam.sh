#!/bin/sh -x
# -*- mode: sh; -*-
# $Id$

# Script to merge existing gene_descriptions with those from the families
# (where known) into a new table called $new_gene_description, which is a
# complete replacement for the old table.  The script assumes that the family
# and ensembl-core database live in the same server (so you can do joins).

usage="Usage: $0 ens_core_database MUSG -h host -u user family_database"
if [ $# -lt 6 ] ; then
    echo $usage
    exit 1
fi

# where mysql lives (or how it can be found, provided PATH is ok):
mysql=mysql; mysql_extra_flags='--batch'
# mysql=cat  # for debugging


# Database to read existing descriptions from:
read_database=$1; shift

# inside the family database, which database to read additional descriptions
# from:
db_name=$1; shift;

# check this:
should_match='^[A-Z][A-Z][A-Z]G$'
if echo "$db_name" | grep -s $should_match > /dev/null 2>&1 ; then
    : # OK
else
    echo "arg 2: database name to use: '$db_name' does not match $should_match" >&2
    exit 2
fi


# Name of the table having the existing descriptions (only changes when the
# schema does, really, but also useful during testing/debugging)
read_gene_desc_table='gene_description'

# the table with new descriptions to be created (in the current database):
new_gene_description='gene_description'

# now produce the SQL (with the $variables  being replaced with their values)
# and pipe this as input into mysql:
(cat <<EOF
SELECT COUNT(*) AS all_old_descriptions
FROM $read_database.$read_gene_desc_table gd
WHERE gd.description IS NOT NULL
  AND gd.description NOT IN ('', 'unknown', 'UNKNOWN');

# creating new local gene_description table;
CREATE TABLE $new_gene_description
SELECT *
FROM  $read_database.$read_gene_desc_table gd 
WHERE gd.gene_id IS NULL;

ALTER TABLE $new_gene_description ADD PRIMARY KEY(gene_id);
# (to avoid duplicates)

# inserting existing desc's:
INSERT INTO $new_gene_description 
  SELECT *
  FROM $read_database.$read_gene_desc_table;

# delete anything that looks remotely unknown:
DELETE FROM $new_gene_description WHERE DESCRIPTION is null;
DELETE FROM $new_gene_description WHERE DESCRIPTION = '';
DELETE FROM $new_gene_description WHERE DESCRIPTION = 'unknown';
DELETE FROM $new_gene_description WHERE DESCRIPTION = 'UNKNOWN';

# selecting all genes from gene, this time properly as 'unknown':
INSERT INTO $new_gene_description
  SELECT id, 'unknown'
  FROM $read_database.gene;
# this will fail on the knowns, as they should; only the 'unknowns' do get
# in.

# give stats:
SELECT COUNT(*) AS unknown_old_descriptions
FROM $new_gene_description
WHERE description = 'unknown';
  
# this table is local to our mysql session, and is deleted automatically
# when it ends
CREATE TEMPORARY TABLE tmp_new_descriptions
  SELECT gd.gene_id, f.description
  FROM family f, family_members fm, $new_gene_description gd
  WHERE gd.description ='unknown'
    AND fm.db_name ='$db_name'
    AND fm.db_id = gd.gene_id
    AND fm.family = f.internal_id
    AND f.description <> 'UNKNOWN';

DELETE FROM $new_gene_description where description = 'unknown';

INSERT INTO $new_gene_description
  SELECT *
  FROM tmp_new_descriptions;

# give new stats:
SELECT COUNT(*) AS all_new_descriptions
FROM $new_gene_description;
EOF
) | $mysql $mysql_extra_flags "$@"
