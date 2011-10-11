#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use DBI qw( :sql_types );
use Getopt::Long qw( :config no_ignore_case );

my $timestamp = strftime( "%Y%m%d-%H%M%S", localtime() );

# Master database location:
my ( $mhost, $mport ) = ( 'ens-staging1', '3306' );
my ( $muser, $mpass ) = ( 'ensro',        undef );
my $mdbname = 'ensembl_production';

# User database location (default values):
my ( $host, $port ) = ( undef, '3306' );
my ( $user, $pass );
my $dbname;
my $dbpattern;

my $core    = 0;
my $verbose = 0;
my $dumppath;

# Do command line parsing.
if ( !GetOptions( 'mhost|mh=s'     => \$mhost,
                  'mport|mP=i'     => \$mport,
                  'muser|mu=s'     => \$muser,
                  'mpass|mp=s'     => \$mpass,
                  'mdatabase|md=s' => \$mdbname,
                  'host|h=s'       => \$host,
                  'port|P=i'       => \$port,
                  'user|u=s'       => \$user,
                  'pass|p=s'       => \$pass,
                  'database|d=s'   => \$dbname,
                  'pattern=s'      => \$dbpattern,
                  'verbose|v!'     => \$verbose,
                  'core=i'         => \$core,
                  'dumppath|dp=s'  => \$dumppath )
     ||
     !( defined($host) &&
        defined($user) &&
        defined($pass) &&
        ( defined($dbname) || defined($dbpattern) || defined($core) ) &&
        defined($mhost) &&
        defined($muser) &&
        defined($dumppath) ) )
{
  my $indent = ' ' x length($0);
  print <<USAGE_END;
This script populates the analysis_description table of a user-defined
database from the production database.

Usage:

  $0 -h host [-P port] \\
  $indent -u user [-p password]
  $indent -d database | --pattern pattern \\
  $indent -dp dumppath \\
  $indent [-s species -t type] \\
  $indent [-mh host] [-mP port] \\
  $indent [-mu user] [-mp password] [-md database] \\
  $indent [-v]

  -h / --host         User database server host.
  -P / --port         User database server port (optional, default is 3306).

  -u / --user         User username (must have write-access).
  -p / --pass         User password.

  -d / --database     User database name or SQL pattern,
                      e.g. --database="homo_sapiens_rnaseq_62_37g"
                      or   --database="%core_62%".

  -dp / --dumppath    Dump path.
                      Back-up table into the specified directory path.

  --pattern           User database by Perl regular expressionm
                      e.g. --pattern="^homo.*(rnaseq|vega)_62".

                      (-d/--database and --pattern are mutually exclusive)

  --core=NN           Preset pattern for Core-like databases in relase NN.
                      Specifying --core=62 is equivalent to using
                      --pattern="(cdna|core|otherfeatures|rnaseq|vega)_62".

  -s / --species      If the name of the database specified by -d (or
                      --database) is not in the standard format, this
                      flag is used to specify what species it is (e.g.,
                      'homo_sapiens', 'gallus_gallus' etc., i.e. the
                      latin name). ** NOT YET IMPLEMENTED **

  -t / --type         If the name of the database specified by -d (or
                      --database) is not in the standard format, this
                      flag may be used to specify what type it is (e.g.,
                      'core', 'vega' etc.). ** NOT YET IMPLEMENTED **

  -mh / --mhost       Production database server host
                      (optional, default is 'ens-staging1').
  -mP / --mport       Production database server port
                      (optional, default is 3306).

  -mu / --muser       Production database username (no write-access required)
                      (optional, default is 'ensro').
  -mp / --mpass       Production database password
                      (optional, default is undefined).

  -md / --mdatabase   Production database name
                      (optional, default is 'ensembl_production').

  -v / --verbose      Be verbose, display every SQL statement as they
                      are executed (on standard error).


USAGE_END

  die( "Need the following options: " .
       "-h -u -p -d (or --pattern) and -dp\n" );

} ## end if ( !GetOptions( 'mhost|mh=s'...))

if ($core) {
  $dbpattern =
    sprintf( '(cdna|core|otherfeatures|rnaseq|vega)_%d', $core );
}

if ( defined($dbname) && defined($dbpattern) ) {
  die("-d/--database and --pattern/--core are mutually exclusive\n");
}

# Fetch all data from the master database.
my %data;
{
  my $dsn = sprintf( 'DBI:mysql:host=%s;port=%d;database=%s',
                     $mhost, $mport, $mdbname );
  my $dbh = DBI->connect( $dsn, $muser, $mpass,
                          { 'PrintError' => 1, 'RaiseError' => 1 } );

  my $sth =
    $dbh->prepare( 'SELECT full_db_name, logic_name, ' .
                  'description, display_label, displayable, web_data ' .
                  'FROM full_analysis_description' );

  $sth->execute();

  my ( $full_db_name, $logic_name, %hash );
  $sth->bind_columns( \( $full_db_name,        $logic_name,
                         $hash{'description'}, $hash{'display_label'},
                         $hash{'displayable'}, $hash{'web_data'} ) );

  while ( $sth->fetch() ) {
    $data{$full_db_name}{$logic_name} = { %{ \%hash } };
  }

  $dbh->disconnect();
}

# Put all data into the specified database.
{
  my $dsn = sprintf( 'DBI:mysql:host=%s;port=%d', $host, $port );
  my $dbh = DBI->connect( $dsn, $user, $pass,
                          { 'PrintError' => 1, 'RaiseError' => 1 } );

  my $sth;
  if ( defined($dbname) ) {
    $sth = $dbh->prepare('SHOW DATABASES LIKE ?');
    $sth->bind_param( 1, $dbname, SQL_VARCHAR );
  }
  else {
    $sth = $dbh->prepare('SHOW DATABASES');
  }

  $sth->execute();

  $sth->bind_col( 1, \$dbname );

  while ( $sth->fetch() ) {
    if ( defined($dbpattern) && $dbname !~ /$dbpattern/ ) { next }

    my $dbdata = $data{$dbname};
    if ( !defined($dbdata) ) {
      printf( "ERROR: Can not find data for database '%s' " .
                "(skipping it)!\n",
              $dbname );
      next;
    }

    print( '=' x 80, "\n" );
    printf( "\t%s\n", $dbname );
    print( '=' x 80, "\n" );

    my $full_table_name =
      $dbh->quote_identifier( undef, $dbname, 'analysis_description' );
    my $full_table_name_bak = $dbh->quote_identifier( undef, $dbname,
                                           'analysis_description_bak' );

    if ( defined($dumppath) ) {
      # Backup the table on file.
      my $filename = sprintf( "%s/%s.%s.%s.sql",
                              $dumppath, $dbname,
                              'analysis_description', $timestamp );

      if ( -e $filename ) {
        die( sprintf( "File '%s' already exists.", $filename ) );
      }

      printf( "Backing up table %s on file.\n",
              'analysis_description' );
      printf( "--> %s\n", $filename );

      if ( system( "mysqldump",
                   "--host=$host",
                   "--port=$port",
                   "--user=$user",
                   ( defined($pass) ? "--password=$pass" : "--opt" ),
                   "--result-file=$filename",
                   "$dbname",
                   'analysis_description' ) )
      {
        die("mysqldump failed: $?");
      }
    } ## end if ( defined($dumppath...))

    $dbh->do(
           sprintf( 'DROP TABLE IF EXISTS %s', $full_table_name_bak ) );

    # Make a backup of any existing data.
    $dbh->do( sprintf( 'CREATE TABLE %s LIKE %s',
                       $full_table_name_bak, $full_table_name ) );
    $dbh->do( sprintf( 'INSERT INTO %s SELECT * FROM %s',
                       $full_table_name_bak, $full_table_name ) );

    # Truncate (empty) the table before inserting new data into it.
    $dbh->do( sprintf( 'TRUNCATE TABLE %s', $full_table_name ) );

    # Get all logic names and their corresponding analysis_id
    my $sth2 = $dbh->prepare(
                sprintf( 'SELECT logic_name, analysis_id FROM %s',
                  $dbh->quote_identifier( undef, $dbname, 'analysis' ) )
    );

    $sth2->execute();

    my ( $logic_name, $analysis_id );
    $sth2->bind_columns( \( $logic_name, $analysis_id ) );

    while ( $sth2->fetch() ) {
      $dbdata->{ lc($logic_name) }{'analysis_id'} = $analysis_id;
    }

    # Insert into the empty analysis_description table.
    $sth2 = $dbh->prepare(
        sprintf( 'INSERT INTO %s ' .
                   '(analysis_id, description, display_label, ' .
                   'displayable, web_data) ' . 'VALUES (?, ?, ?, ?, ?)',
                 $full_table_name ) );

    foreach my $logic_name ( keys( %{$dbdata} ) ) {
      if ( !exists( $dbdata->{$logic_name}{'description'} ) ) {
        printf( "ERROR: Missing production database entry " .
                  "for logic name '%s'\n",
                $logic_name );
      }
      elsif ( !exists( $dbdata->{$logic_name}{'analysis_id'} ) ) {
        printf( "WARNING: Expected to find analysis entry " .
                  "for logic name '%s'\n",
                $logic_name );
      }
      else {
        $sth2->bind_param( 1, $dbdata->{$logic_name}{'analysis_id'},
                           SQL_INTEGER );
        $sth2->bind_param( 2, $dbdata->{$logic_name}{'description'},
                           SQL_VARCHAR );
        $sth2->bind_param( 3, $dbdata->{$logic_name}{'display_label'},
                           SQL_VARCHAR );
        $sth2->bind_param( 4, $dbdata->{$logic_name}{'displayable'},
                           SQL_INTEGER );
        $sth2->bind_param( 5, $dbdata->{$logic_name}{'web_data'},
                           SQL_VARCHAR );

        $sth2->execute();
      }
    } ## end foreach my $logic_name ( keys...)

    print("Inserted data\n");

    my $key_name = 'analysis_id';
    {
      my $statement = sprintf( 'SELECT %s ' . 'FROM %s ' .
                                 'LEFT JOIN %s t USING (%s) ' .
                                 'WHERE t.%s IS NULL ' . 'ORDER BY %s',
                               $key_name,            $full_table_name,
                               $full_table_name_bak, $key_name,
                               $key_name,            $key_name );

      my $sth3 = $dbh->prepare($statement);

      $sth3->execute();

      my $key;
      $sth3->bind_col( 1, \$key );

      my @keys;
      while ( $sth3->fetch() ) {
        push( @keys, $key );
      }

      if (@keys) {
        print("\n");
        print("New data inserted:\n");
        printf( "SELECT * FROM %s WHERE %s IN (%s);\n",
                $full_table_name, $key_name, join( ',', @keys ) );
        print("\n");
      }
    }
    {
      my $statement = sprintf( 'SELECT %s ' . 'FROM %s ' .
                                 'LEFT JOIN %s t USING (%s) ' .
                                 'WHERE t.%s IS NULL ' . 'ORDER BY %s',
                               $key_name,        $full_table_name_bak,
                               $full_table_name, $key_name,
                               $key_name,        $key_name );

      my $sth3 = $dbh->prepare($statement);

      $sth3->execute();

      my $key;
      $sth3->bind_col( 1, \$key );

      my @keys;
      while ( $sth3->fetch() ) {
        push( @keys, $key );
      }

      if (@keys) {
        print("\n");
        print( '-' x 40, "\n" );
        print("!! Old data deleted:\n");
        printf( "SELECT * FROM %s WHERE %s IN (%s);\n",
                $full_table_name_bak, $key_name, join( ',', @keys ) );
        print( '-' x 40, "\n" );
        print("\n");
      }

    }
    #delete the backup table
    $dbh->do(
           sprintf( 'DROP TABLE IF EXISTS %s', $full_table_name_bak ) );

  } ## end while ( $sth->fetch() )
  continue {
    print("\n");
  }

  $dbh->disconnect();
}

print <<FINAL_END
To restore a table from dump login to the database and use command:
"source {dumpfile} or "\. {dumpfile}".
FINAL_END
