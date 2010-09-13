#!/usr/bin/env perl

use strict;
use warnings;

use DBI qw(:sql_types);

my @staging_servers = ( { 'dbhost' => 'ens-staging1',
                          'dbport' => '3306' }, {
                          'dbhost' => 'ens-staging2',
                          'dbport' => '3306' } );

my $master_server      = $staging_servers[0];
my $production_db_name = 'ensembl_production_60';

my $dbuser   = 'ensadmin';
my $dbpass   = 'ensembl';
my $dbrouser = 'ensro';
my $dbropass = undef;

my %databases;
foreach my $server (@staging_servers) {
  my $dsn = sprintf( 'DBI:mysql:host=%s;port=%s',
                     $server->{'dbhost'}, $server->{'dbport'} );
  my $dbh = DBI->connect( $dsn, $dbrouser, $dbropass,
                          { 'PrintError' => 1, 'RaiseError' => 0 } );

  foreach my $dbtype ( 'cdna',                'core',
                       'coreexpressionatlas', 'coreexpressionest',
                       'coreexpressiongnf',   'funcgen',
                       'otherfeatures',       'variation',
                       'vega' )
  {
    my $sth = $dbh->prepare(
             sprintf( "SHOW DATABASES LIKE '%%\\_%s\\_%%'", $dbtype ) );

    $sth->execute();

    my $dbname;
    $sth->bind_col( 1, \$dbname );

    while ( $sth->fetch() ) {
      if ( $dbname !~
           /^([a-z]+_[a-z]+)_([a-z]+)_([0-9]+)_([0-9]+)([a-z]?)$/
           || exists( $databases{$1}{$2} ) )
      {
        next;
      }

      if ( $2 ne $dbtype ) {
        carp( sprintf( "Strange database type '%s', expected '%s'",
                       $2, $dbtype ) );
        next;
      }

      my $cn_sth = $dbh->prepare(
                         sprintf(
                           "SELECT meta_value "
                             . "FROM %s.meta "
                             . "WHERE meta_key = 'species.common_name'",
                           $dbh->quote_identifier($dbname) ) );

      $cn_sth->execute();

      my $common_name;
      if ( !$cn_sth->err() ) {
        $cn_sth->bind_col( 1, \$common_name );
        while ( $cn_sth->fetch() ) { }
      }

      $databases{$1}{$2} = { 'db_release'  => $3,
                             'db_assembly' => $4,
                             'db_suffix'   => $5,
                             'db_host'     => $server->{'dbhost'},
                             'common_name' => $common_name };

    } ## end while ( $sth->fetch() )
  } ## end foreach my $dbtype ( 'cdna'...)

  $dbh->disconnect();

} ## end foreach my $server (@staging_servers)

my $dsn = sprintf( 'DBI:mysql:host=%s;port=%s;database=%s',
                   $master_server->{'dbhost'},
                   $master_server->{'dbport'},
                   $production_db_name );
my $dbh = DBI->connect( $dsn, $dbuser, $dbpass,
                        { 'PrintError' => 0, 'RaiseError' => 0 } );

my $sp_sel_sth = $dbh->prepare(
  q(
SELECT species_id
FROM species
WHERE db_name = ?
) );

my $sp_sth = $dbh->prepare(
  q(
INSERT INTO species
  (db_name, common_name)
VALUES (?, ?)
) );

my $db_sth = $dbh->prepare(
  q(
INSERT INTO db
  (species_id, db_type, db_release, db_assembly, db_suffix, db_host)
VALUES  (?, ?, ?, ?, ?, ?)
) );

foreach my $db_name ( sort( keys(%databases) ) ) {
  foreach my $db_type ( sort ( keys( %{ $databases{$db_name} } ) ) ) {
    my $entry = $databases{$db_name}{$db_type};

    $sp_sel_sth->bind_param( 1, $db_name, SQL_VARCHAR );

    $sp_sel_sth->execute();

    my $species_id;
    $sp_sel_sth->bind_col( 1, \$species_id );

    while ( $sp_sel_sth->fetch() ) { }

    if ( !defined($species_id) ) {
      $sp_sth->bind_param( 1, $db_name, SQL_VARCHAR );
      $sp_sth->bind_param( 2, $entry->{'common_name'}, SQL_VARCHAR );

      printf( "Inserting '%s' ('%s') into species_list table... ",
              $db_name, $entry->{'common_name'} );

      $sp_sth->execute();

      if ( $sp_sth->err() ) {
        print("failed\n");
        next;
      } else {
        print("ok\n");
      }

      $species_id = $dbh->{'mysql_insertid'};
    }

    if ( !defined($species_id) ) {
      die( sprintf( "No species_id for '%s'.", $db_name ) );
    }

    $db_sth->bind_param( 1, $species_id,             SQL_INTEGER );
    $db_sth->bind_param( 2, $db_type,                SQL_VARCHAR );
    $db_sth->bind_param( 3, $entry->{'db_release'},  SQL_INTEGER );
    $db_sth->bind_param( 4, $entry->{'db_assembly'}, SQL_INTEGER );
    $db_sth->bind_param( 5, $entry->{'db_suffix'},   SQL_VARCHAR );
    $db_sth->bind_param( 6, $entry->{'db_host'},     SQL_VARCHAR );

    printf( "Inserting database '%s_{%s}_%d_%d%s'... ",
            $db_name, $db_type,
            $entry->{'db_release'},
            $entry->{'db_assembly'},
            $entry->{'db_suffix'} );

    $db_sth->execute();

    if ( $db_sth->err() ) {
      print("failed\n");
    } else {
      print("ok\n");
    }
  } ## end foreach my $db_type ( sort ...)

} ## end foreach my $db_name ( sort(...))

$dbh->disconnect();
