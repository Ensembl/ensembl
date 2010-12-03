#!/usr/bin/env perl

# This little script will bootstrap the master tables in the production
# database (and nothing else).  This means copying the relevant tables
# from a known correct database (the "template database") into the
# master_% tables in the production database.
#
# The template database needs to live on the same server as the
# production database.
#
# The output of this script is SQL written to standard output.  The SQL
# needs to be run against the MySQL server manually.
#

use strict;
use warnings;

my $template_db = $ARGV[0];

if ( !defined($template_db) ) {
  print STDERR <<END_USAGE;
Usage:
  $0 template_database_name >output_file.sql

  "template_database_name" should be the name of a Core database.
END_USAGE
  exit;
}

#-----------------------------------------------------------------------
# The "simple" tables.
{

  my @simple_tables =
    ( 'attrib_type', 'external_db', 'misc_set', 'unmapped_reason' );

  foreach my $table (@simple_tables) {
    print(
      qq(
DROP TABLE IF EXISTS master_${table};
CREATE TABLE master_${table}
  LIKE  ${template_db}.${table};

INSERT INTO master_${table}
  SELECT  *
    FROM  ${template_db}.${table};
) );
  }

}
