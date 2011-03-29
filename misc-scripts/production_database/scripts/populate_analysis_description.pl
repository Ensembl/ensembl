#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use DBI qw( :sql_types );
use Getopt::Long qw( :config no_ignore_case );

local $| = 1;

# Master database location:
my ( $mhost, $mport ) = ( 'ens-staging1', '3306' );
my ( $muser, $mpass ) = ( 'ensro',        undef );
my $mdbname = 'ensembl_production';

# User database location (default values):
my ( $host, $port ) = ( undef, '3306' );
my ( $user, $pass );
my $dbname;
my $dbpattern;

my $core     = 0;
my $drop_bak = 0;
my $verbose  = 0;
my $cleanup  = 0;

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
                  'drop|D!'        => \$drop_bak,
                  'verbose|v!'     => \$verbose,
                  'core=i'         => \$core,
                  'cleanup!'       => \$cleanup )
     || !(
           defined($host)
        && defined($user)
        && defined($pass)
        && ( defined($dbname) || defined($dbpattern) || defined($core) )
        && defined($mhost)
        && defined($muser) ) )
{
  my $indent = ' ' x length($0);
  print <<USAGE_END;
This script populates the analysis_description table of a user-defined
database from the production database.

Usage:

  $0 -h host [-P port] \\
  $indent -u user [-p password]
  $indent -d database | --pattern pattern \\
  $indent [-mh host] [-mP port] \\
  $indent [-mu user] [-mp password] [-md database] \\
  $indent [-D] [-v]

  -h / --host       User database server host
  -P / --port       User database server port (optional, default is 3306)

  -u / --user       User username (must have write-access)
  -p / --pass       User password

  -d / --database   User database name or SQL pattern
                    e.g. --database="homo_sapiens_rnaseq_62_37g"
                    or   --database="%core_62%"

  --pattern         User database by Perl regular expression
                    e.g. --pattern="^homo.*(rnaseq|vega)_62"

                    (-d/--database and --pattern are mutually exclusive)

  --core=NN         Preset pattern for Core-like databases in relase NN
                    Specifying --core=62 is equivalent to using
                    --pattern="(cdna|core|otherfeatures|rnaseq|vega)_62"

  -mh / --mhost     Production database server host
                    (optional, default is 'ens-staging1')
  -mP / --mport     Production database server port
                    (optional, default is 3306)

  -mu / --muser     Production database username (no write-access required)
                    (optional, default is 'ensro')
  -mp / --mpass     Production database password
                    (optional, default is undefined)

  -md / --mdatabase Production database name
                    (optional, default is 'ensembl_production')

  -D / --drop       Drop backup tables if they already exists in the
                    database from a previous run

  --cleanup         Just clean up backup tables, do nothing else

  -v / --verbose    Be verbose, display every SQL statement as they
                    are executed (on standard error)

USAGE_END

  die(   "Need the following options: "
       . "-h -u -p and -d (or --pattern)\n" );

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
    $dbh->prepare(  'SELECT full_db_name, logic_name, '
                  . 'description, display_label, displayable, web_data '
                  . 'FROM full_analysis_description' );

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
  } else {
    $sth = $dbh->prepare('SHOW DATABASES');
  }

  $sth->execute();

  $sth->bind_col( 1, \$dbname );

  while ( $sth->fetch() ) {
    if ( defined($dbpattern) && $dbname !~ /$dbpattern/ ) { next }

    my $dbdata = $data{$dbname};
    if ( !defined($dbdata) ) {
      printf( "ERROR: Can not find data for database '%s' "
                . "(skipping it)!\n",
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

    # Drop backup table if it exists and if asked to do so.
    if ( $cleanup || $drop_bak ) {
      $dbh->do(
           sprintf( 'DROP TABLE IF EXISTS %s', $full_table_name_bak ) );

      print("Dropped old backup table\n");
      if ($cleanup) { next }
    }

    # Make a backup of any existing data.
    $dbh->do( sprintf( 'CREATE TABLE %s LIKE %s',
                       $full_table_name_bak, $full_table_name ) );
    $dbh->do( sprintf( 'INSERT INTO %s SELECT * FROM %s',
                       $full_table_name_bak, $full_table_name ) );
    print("Created new backup table\n");

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
               sprintf( 'INSERT INTO %s '
                          . '(analysis_id, description, display_label, '
                          . 'displayable, web_data) '
                          . 'VALUES (?, ?, ?, ?, ?)',
                        $full_table_name ) );

    foreach my $logic_name ( keys( %{$dbdata} ) ) {
      if ( !exists( $dbdata->{$logic_name}{'description'} ) ) {
        printf( "ERROR: Missing production database entry "
                  . "for logic name '%s'\n",
                $logic_name );
      } elsif ( !exists( $dbdata->{$logic_name}{'analysis_id'} ) ) {
        printf( "WARNING: Expected to find analysis entry "
                  . "for logic name '%s'\n",
                $logic_name );
      } else {
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
    }

    print("Inserted data\n");

    my $key_name = 'analysis_id';
    {
      my $statement = sprintf( 'SELECT %s '
                                 . 'FROM %s '
                                 . 'LEFT JOIN %s t USING (%s) '
                                 . 'WHERE t.%s IS NULL '
                                 . 'ORDER BY %s',
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
      my $statement = sprintf( 'SELECT %s '
                                 . 'FROM %s '
                                 . 'LEFT JOIN %s t USING (%s) '
                                 . 'WHERE t.%s IS NULL '
                                 . 'ORDER BY %s',
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

      print("To undo:\n");
      printf( "  DROP TABLE %s;\n", $full_table_name );
      printf( "  RENAME TABLE %s TO %s;\n",
              $full_table_name_bak, $full_table_name );

    }
  } continue {
    print("\n");
  }

  if ( !$cleanup ) {
    print("Remember to drop the backup tables!\n\n");
  }

  $dbh->disconnect();
}
