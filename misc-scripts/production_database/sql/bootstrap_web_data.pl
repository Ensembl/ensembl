#!/usr/bin/env perl

# This little script will bootstrap the web_data table in the production
# database (and nothing else).
#
# This will be done by using the db table (populated by the
# update_database_list.pl script) to iterate over current databases.
#
# The output of this script is SQL written directly to the production
# database.  Use with caution.

use strict;
use warnings;

use DBI (':sql_types');
use Data::Dumper;

$Data::Dumper::Terse    = 1;
$Data::Dumper::Useqq    = 0;
$Data::Dumper::Indent   = 0;
$Data::Dumper::Deparse  = 0;
$Data::Dumper::Sortkeys = 1;

my $dbuser = 'XXX';
my $dbpass = 'YYY';

my $dsn = sprintf( "DBI:mysql:host=%s;port=%d;database=%s",
                   'ens-staging1', 3306, 'ensembl_production',
                   { 'PrintError' => 1 } );
my $dbh = DBI->connect( $dsn, $dbuser, $dbpass );

my %servers = (
      'ens-staging1' => {
        'dbh' =>
          DBI->connect(
          sprintf( "DBI:mysql:host=%s;port=%d;", 'ens-staging1', 3306 ),
          $dbuser, $dbpass
          )
      },
      'ens-staging2' => {
        'dbh' =>
          DBI->connect(
          sprintf( "DBI:mysql:host=%s;port=%d;", 'ens-staging2', 3306 ),
          $dbuser, $dbpass
          ) } );

my $sth = $dbh->prepare(
  q(
SELECT db_list.full_db_name, db.db_host, db.species_id, db.db_type
FROM db_list
  JOIN db USING (db_id)
WHERE db.is_current = 1
AND db.db_type IN ( 'core', 'otherfeatures', 'cdna', 'vega', 'rnaseq' )
) );

$sth->execute();

my ( $db_name, $db_host, $species_id, $db_type );
$sth->bind_columns( \( $db_name, $db_host, $species_id, $db_type ) );

my $sth3 = $dbh->prepare(
  q(
SELECT analysis_description_id
FROM analysis_description
WHERE logic_name = ?
) );

my $sth4 = $dbh->prepare(
  q(
INSERT INTO web_data (web_data_id, data) VALUE (?,?)
) );

my $sth5 = $dbh->prepare(
  q(
INSERT INTO analysis_web_data (
  analysis_description_id,
  web_data_id,
  species_id,
  db_type,
  displayable
) VALUES (?, ?, ?, ?, ?)
) );

my %web_data;

while ( $sth->fetch() ) {
  printf( "%s@%s\n", $db_name, $db_host );
  my $sth2 = $servers{$db_host}{'dbh'}->prepare(
    sprintf(
      q(
SELECT a.logic_name, ad.web_data, ad.displayable
FROM %s a
  JOIN %s ad USING (analysis_id)
),
      $dbh->quote_identifier( undef, $db_name, 'analysis' ),
      $dbh->quote_identifier( undef, $db_name, 'analysis_description'
      ),
    ) );

  $sth2->execute();

  my ( $logic_name, $web_data, $displayable );
  $sth2->bind_columns( \( $logic_name, $web_data, $displayable ) );

  while ( $sth2->fetch() ) {

    $sth3->bind_param( 1, $logic_name, SQL_VARCHAR );
    $sth3->execute();

    my $id;
    $sth3->bind_col( 1, \$id );

    while ( $sth3->fetch() ) {
      my $key =
        ( defined($web_data) ? Dumper( eval($web_data) ) : 'NULL' );
      push( @{ $web_data{$key} }, {
              'analysis_description_id' => $id,
              'species_id'              => $species_id,
              'db_type'                 => $db_type,
              'displayable'             => $displayable
            } );
      last;
    }

    $sth3->finish();

  }
} ## end while ( $sth->fetch() )

my $id = 0;
foreach my $web_data ( keys(%web_data) ) {
  if ( $web_data ne 'NULL' ) {
    $sth4->bind_param( 1, ++$id,     SQL_INTEGER );
    $sth4->bind_param( 2, $web_data, SQL_VARCHAR );
    $sth4->execute();
  }

  foreach my $data ( @{ $web_data{$web_data} } ) {
    $sth5->bind_param( 1, $data->{'analysis_description_id'},
                       SQL_INTEGER );
    if ( $web_data eq 'NULL' ) {
      $sth5->bind_param( 2, undef, SQL_INTEGER );
    } else {
      $sth5->bind_param( 2, $id, SQL_INTEGER );
    }
    $sth5->bind_param( 3, $data->{'species_id'},  SQL_INTEGER );
    $sth5->bind_param( 4, $data->{'db_type'},     SQL_VARCHAR );
    $sth5->bind_param( 5, $data->{'displayable'}, SQL_INTEGER );
    $sth5->execute();
  }
}

