#!/usr/bin/perl -w

use strict;
use warnings;

use DBI qw( :sql_types );
use Getopt::Long qw( :config no_ignore_case );

#-----------------------------------------------------------------------

sub usage {
  print("Usage:\n");
  printf( "\t%s\t-h dbhost [-P dbport] \\\n"
      . "\t%s\t-u dbuser [-p dbpass] \\\n"
      . "\t%2\$s\t-d dbname\n",
    $0, ' ' x length($0) );
  print("\n");
  printf( "\t%s\t-?\n", $0 );
  print("\n");
  print("Arguments:\n");
  print("\t-h/--host dbhost\tDatabase server host name\n");
  print("\t-P/--port dbport\tDatabase server port (optional)\n");
  print("\t-u/--user dbuser\tDatabase user name\n");
  print("\t-p/--pass dbpass\tUser password (optional)\n");
  print("\t-d/--name dbname\tDatabase name\n");
  print("\t-?/--help\t\tDisplays this information\n");
}

#-----------------------------------------------------------------------

my ( $dbhost, $dbport );
my ( $dbuser, $dbpass );
my $dbname;

$dbport = '3306';

if (
  !GetOptions(
    'dbhost|host|h=s' => \$dbhost,
    'dbport|port|P=i' => \$dbport,
    'dbuser|user|u=s' => \$dbuser,
    'dbpass|pass|p=s' => \$dbpass,
    'dbname|name|d=s' => \$dbname,
    'help|?'          => sub { usage(); exit } )
  || !defined($dbhost)
  || !defined($dbuser)
  || !defined($dbname) )
{
  usage();
  exit;
}

my $dsn = sprintf( "DBI:mysql:database=%s;host=%s;port=%s",
  $dbname, $dbhost, $dbport );

my $dbh =
  DBI->connect( $dsn, $dbuser, $dbpass,
  { 'RaiseError' => 0, 'PrintError' => 0 } );

my %subsets;
{
  my $statement = q(
SELECT DISTINCT
        ontology.name,
        subset.name
FROM    ontology,
        term,
        subset
WHERE   ontology.ontology_id = term.ontology_id
  AND   FIND_IN_SET(subset.name, term.subsets) > 0
);

  my $sth = $dbh->prepare($statement);

  $sth->execute();

  my ( $ontology_name, $subset_name );

  $sth->bind_columns( \( $ontology_name, $subset_name ) );

  while ( $sth->fetch() ) {
    push( @{ $subsets{$ontology_name} }, $subset_name );
  }

  $sth->finish();
}

{
  my $select_statement = q(
SELECT  child_term.term_id,
        parent_term.term_id,
        closure.distance
FROM    closure,
        term child_term,
        term parent_term,
        ontology
WHERE   closure.parent_term_id = parent_term.term_id
  AND   closure.child_term_id = child_term.term_id
  AND   FIND_IN_SET(?, parent_term.subsets) > 0
  AND   parent_term.ontology_id = ontology.ontology_id
  AND   ontology.name = ?
ORDER BY child_term.accession, closure.distance;
);

  my $select_sth = $dbh->prepare($select_statement);

  foreach my $ontology_name ( keys(%subsets) ) {
    foreach my $subset_name ( @{ $subsets{$ontology_name} } ) {

      my $aux_table_name = $dbh->quote_identifier(
        sprintf( "aux_%s_%s_map", $ontology_name, $subset_name ) );

      $dbh->do(
        sprintf(
          "CREATE TABLE %s ( "
            . "term_id INT UNSIGNED NOT NULL, "
            . "subset_term_id INT UNSIGNED NOT NULL, "
            . "UNIQUE INDEX map_idx (term_id, subset_term_id) )",
          $aux_table_name
        ) );

      if ( $dbh->err() ) {
        printf( "MySQL error, \"%s\", skipping...\n", $dbh->errstr() );
        next;
      }

      $select_sth->bind_param( 1, $subset_name,   SQL_VARCHAR );
      $select_sth->bind_param( 2, $ontology_name, SQL_VARCHAR );

      $select_sth->execute();

      my ( $child_id, $parent_id, $distance );

      $select_sth->bind_columns(
        \( $child_id, $parent_id, $distance ) );

      my $insert_statement = sprintf(
        "INSERT IGNORE INTO %s "
          . "(term_id, subset_term_id) "
          . "VALUES (?, ?)",
        $aux_table_name
      );

      $dbh->do( sprintf( "LOCK TABLE %s WRITE", $aux_table_name ) );

      my $insert_sth = $dbh->prepare($insert_statement);

      printf( "%s...\n", $aux_table_name );

      my $last_child_id;
      my $the_distance;

      while ( $select_sth->fetch() ) {
        if ( !defined($last_child_id)
          || $child_id != $last_child_id )
        {
          $last_child_id = $child_id;
          $the_distance  = $distance;
        }

        if ( $child_id == $last_child_id
          && $distance != $the_distance )
        {
          next;
        }

        if ( $child_id == $last_child_id
          && $distance == $the_distance )
        {
          $insert_sth->bind_param( 1, $child_id,  SQL_INTEGER );
          $insert_sth->bind_param( 2, $parent_id, SQL_INTEGER );

          $insert_sth->execute();
        }
      }

      $select_sth->finish();

      $dbh->do( sprintf( "OPTIMIZE TABLE %s", $aux_table_name ) );
      $dbh->do("UNLOCK TABLES");

    } ## end foreach my $subset_name ( @...)
  } ## end foreach my $ontology_name (...)
}

# $Id$
