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
my @servers = ( 'ens-staging1', 'ens-staging2' );
my $master = 'ens-staging1';

my $dbport = '3306';
my ( $dbwuser, $dbwpass ) = ( 'ensadmin', undef );
my ( $dbuser,  $dbpass )  = ( 'ensro',    undef );

my $opt_help  = 0;
my $opt_about = 0;

if ( !GetOptions( 'release|r=i'  => \$release,
                  'master|m=s'   => \$master,
                  'server|s=s@'  => \@servers,
                  'dbuser|u=s'   => \$dbuser,
                  'dbpass|p=s'   => \$dbpass,
                  'dbport|P=s'   => \$dbport,
                  'dbwuser|wu=s' => \$dbwuser,
                  'dbwpass|wp=s' => \$dbwpass,
                  'help|h!'      => \$opt_help,
                  'about!'       => \$opt_about )
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

my %species;

{
  my $dsn = sprintf( 'DBI:mysql:host=%s;port=%d;database=%s',
                     $master, $dbport, 'ensembl_production' );
  my $dbh = DBI->connect( $dsn, $dbuser, $dbpass,
                          { 'PrintError' => 1, 'RaiseError' => 1 } );

  my $statement =
      'SELECT species_id, web_name, db_name '
    . 'FROM species '
    . 'WHERE is_current = 1';

  my $sth = $dbh->prepare($statement);
  $sth->execute();

  my ( $species_id, $web_name, $db_name );
  $sth->bind_columns( \( $species_id, $web_name, $db_name ) );

  while ( $sth->fetch() ) {
    $species{$db_name} = { 'species_id' => $species_id,
                           'db_name'    => $db_name,
                           'web_name'   => $web_name };
  }

  $dbh->disconnect();
}

my %databases;

foreach my $server (@servers) {
  my $dsn = sprintf( 'DBI:mysql:host=%s;port=%d', $server, $dbport );
  my $dbh = DBI->connect( $dsn, $dbuser, $dbpass,
                          { 'PrintError' => 1, 'RaiseError' => 0 } );

  my $statement = 'SHOW DATABASES LIKE ?';

  my $sth = $dbh->prepare($statement);

  foreach my $species ( keys(%species) ) {
    $sth->bind_param( 1, sprintf( '%s%%\_%d\_%%', $species, $release ),
                      SQL_VARCHAR );

    $sth->execute();

    my $database;

    $sth->bind_col( 1, \$database );

    while ( $sth->fetch() ) {
      my ( $db_type, $db_assembly, $db_suffix ) = ( $database =~
                   /^[a-z]+_[a-z]+_([a-z]+)_[0-9]+_([0-9]+)([a-z]?)$/ );

      if (    !defined($db_type)
           || !defined($db_assembly)
           || !defined($db_suffix) )
      {
        die(
           sprintf( "Failed to parse database name '%s'", $database ) );
      }

      $databases{$database} = {
                       'species_id' => $species{$species}{'species_id'},
                       'db_type'    => $db_type,
                       'db_assembly' => $db_assembly,
                       'db_suffix'   => $db_suffix,
                       'db_host'     => $server };
    }
  } ## end foreach my $species ( keys(...))

  $dbh->disconnect();
} ## end foreach my $server (@servers)

{
  my $dsn = sprintf( 'DBI:mysql:host=%s;port=%d;database=%s',
                     $master, $dbport, 'ensembl_production' );
  my $dbh = DBI->connect( $dsn, $dbwuser, $dbwpass,
                          { 'PrintError' => 1, 'RaiseError' => 1 } );

  my $select_stmt = 'SELECT db_id FROM db_list WHERE full_db_name = ?';
  my $select_sth  = $dbh->prepare($select_stmt);

  my $insert_stmt =
      'INSERT INTO db '
    . '(species_id, db_type, db_release, db_assembly, db_suffix, db_host) '
    . 'VALUES (?, ?, ?, ?, ?, ?)';
  my $insert_sth = $dbh->prepare($insert_stmt);

  foreach my $database ( keys(%databases) ) {
    $select_sth->bind_param( 1, $database, SQL_VARCHAR );
    $select_sth->execute();

    my $db_id;
    $select_sth->bind_col( 1, \$db_id );

    my $is_found = 0;
    while ( $select_sth->fetch() ) {
      $is_found = 1;
      last;
    }
    $select_sth->finish();

    if ($is_found) {
      # printf( "The database '%s' is already in "
      #           . "the production database (id %d)\n",
      #         $database, $db_id );
    } else {
      my $db = $databases{$database};

      printf( "Inserting database '%s' into "
                . "the production database\n",
              $database );

      $insert_sth->bind_param( 1, $db->{'species_id'},  SQL_INTEGER );
      $insert_sth->bind_param( 2, $db->{'db_type'},     SQL_VARCHAR );
      $insert_sth->bind_param( 3, $release,             SQL_INTEGER );
      $insert_sth->bind_param( 4, $db->{'db_assembly'}, SQL_INTEGER );
      $insert_sth->bind_param( 5, $db->{'db_suffix'},   SQL_VARCHAR );
      $insert_sth->bind_param( 6, $db->{'db_host'},     SQL_VARCHAR );

      $insert_sth->execute();
    }
  } ## end foreach my $database ( keys...)

  $dbh->disconnect();
}
