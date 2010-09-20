#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Long qw( :config no_ignore_case );
use DBI qw( :sql_types );

sub usage {
  my $padding = ' ' x length($0);

  print <<USAGE_END;
Usage:
  $0 --release NN --master master-server \\
  $padding --server server1 --server server2 [...] \\
  $padding --dbport 3306 --dbuser user --dbpass passwd \\
  $padding --dbwuser write_user --dbwpass write_passwd

or
  $0 --help

or
  $0 --about

where

  --release/-r  The current release (required).

  --master/-m   The master server where the production database lives
                (optional, default is 'ens-staging1').
  --server/-s   A database server (optional, may occur several times,
                default is 'ens-staging1', and 'ens-staging2').

  --dbport/-P   The port to connect to (optional, default is '3306').

  --dbuser/-u   The (read only) user to connect as (optional,
                default is 'ensro').
  --dbpass/-p   The password to connect with as the above user
                (optional, no default).

  --dbwuser/-wu The user (with write permissions) to connect as
                (optional, default is 'ensadmin').
  --dbwpass/-wp The password to connect with as the above user
                (optional, no default).

  --help/-h     Displays this help text.

  --about/-a    Displays a text about this program (what it does etc.).
USAGE_END
} ## end sub usage

my $release;
my @servers = ('ens-staging1', 'ens-staging2');
my $master = 'ens-staging1';

my $dbport = '3306';
my $dbwuser = 'ensadmin';
my $dbwpass;
my $dbuser = 'ensro';
my $dbpass;

my $opt_help  = 0;
my $opt_about = 0;

if ( !GetOptions( 'release|r=i' => \$release,
                  'master|m=s'  => \$master,
                  'server|s=s@' => \@servers,
                  'dbuser|u=s'  => \$dbuser,
                  'dbpass|p=s'  => \$dbpass,
                  'dbport|P=s'  => \$dbport,
                  'dbrouser|wu' => \$dbwuser,
                  'dbropass|wp' => \$dbwpass,
                  'help|h!'     => \$opt_help,
                  'about!'      => \$opt_about )
     || $opt_help )
{
  usage();
  exit();
} elsif ($opt_about) {
  about();
  exit();
} elsif ( !defined($release) ) {
  print("ERROR: Release was not specified! (use -r or --release)\n");
  usage();
  exit();
}

my %databases;
foreach my $server (@servers) {
  my $dsn = sprintf( 'DBI:mysql:host=%s;port=%d', $server, $dbport );
  my $dbh = DBI->connect( $dsn, $dbuser, $dbpass,
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
                             'db_host'     => $server,
                             'common_name' => $common_name };

    } ## end while ( $sth->fetch() )
  } ## end foreach my $dbtype ( 'cdna'...)

  $dbh->disconnect();

} ## end foreach my $server (@servers)

die;

my $dsn = sprintf( 'DBI:mysql:host=%s;port=%s;database=%s',
                   $master, $dbport,
                   sprintf( 'ensembl_production_%d', $release ) );

my $dbh = DBI->connect( $dsn, $dbwuser, $dbwpass,
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
