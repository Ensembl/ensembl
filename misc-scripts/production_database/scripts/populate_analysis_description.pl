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
my $verbose  = 0;
my $dumppath;
my $dbname_file;

# Do command line parsing.
if ( !GetOptions( 'mhost|mh=s'        => \$mhost,
                  'mport|mP=i'        => \$mport,
                  'muser|mu=s'        => \$muser,
                  'mpass|mp=s'        => \$mpass,
                  'mdatabase|md=s'    => \$mdbname,
                  'host|h=s'          => \$host,
                  'port|P=i'          => \$port,
                  'user|u=s'          => \$user,
                  'pass|p=s'          => \$pass,
                  'database|d=s'      => \$dbname,
                  'pattern=s'         => \$dbpattern,
                  'verbose|v!'        => \$verbose,
                  'core=i'            => \$core,
                  'dumppath|dp=s'     => \$dumppath,
                  'dbname_file|dbf=s' => \$dbname_file )
     || !(
           defined($host)
        && defined($user)
        && defined($pass)
        && ( defined($dbname) || defined($dbpattern) || defined($core) )
        && defined($mhost)
        && defined($muser)
        && defined($dumppath) ) )
{
  my $indent = ' ' x length($0);
  print <<USAGE_END;
This script populates the analysis_description table of a user-defined
database from the production database.

Usage:

  $0 -h host [-P port] \\
  $indent -u user [-p password]
  $indent -d database | --pattern pattern \\
  $indent -dp dumppath
  $indent [-dbf dbname override file]
  $indent [-mh host] [-mP port] \\
  $indent [-mu user] [-mp password] [-md database] \\
  $indent [-v]

  -h / --host         User database server host
  -P / --port         User database server port (optional, default is 3306)

  -u / --user         User username (must have write-access)
  -p / --pass         User password

  -d / --database     User database name or SQL pattern
                      e.g. --database="homo_sapiens_rnaseq_62_37g"
                      or   --database="%core_62%"

  -dp / --dumppath    Dump path.
                      Dump out table into the specified directory path.

  --pattern           User database by Perl regular expression
                      e.g. --pattern="^homo.*(rnaseq|vega)_62"

                      (-d/--database and --pattern are mutually exclusive)

  --core=NN           Preset pattern for Core-like databases in relase NN
                      Specifying --core=62 is equivalent to using
                      --pattern="(cdna|core|otherfeatures|rnaseq|vega)_62"

  -dbf / --dbname_file  Override database names stored in the production
                        database with database names in the specified
                        file.  The file should have the following
                        tab/space delimited data on each line:

                        - database name stored in 'full_db_name' column
                          'db_list' table in the production database.

                        - new database name to use.

  -mh / --mhost       Production database server host
                      (optional, default is 'ens-staging1')
  -mP / --mport       Production database server port
                      (optional, default is 3306)

  -mu / --muser       Production database username (no write-access required)
                      (optional, default is 'ensro')
  -mp / --mpass       Production database password
                      (optional, default is undefined)

  -md / --mdatabase   Production database name
                      (optional, default is 'ensembl_production')

  -v / --verbose      Be verbose, display every SQL statement as they
                      are executed (on standard error)

USAGE_END

  die(   "Need the following options: "
       . "-h -u -p -d (or --pattern) and -dp\n" );

} ## end if ( !GetOptions( 'mhost|mh=s'...))

if ($core) {
  $dbpattern =
    sprintf( '(cdna|core|otherfeatures|rnaseq|vega)_%d', $core );
}

if ( defined($dbname) && defined($dbpattern) ) {
  die("-d/--database and --pattern/--core are mutually exclusive\n");
}

my %dbname_override;

if ( defined($dbname_file) ) {
  if ( !-r $dbname_file ) {
    die(
       sprintf( "File '%s' not found or not readable", $dbname_file ) );
  }

  open( DBNAMEF, "<$dbname_file" );

  my @temp_array;
  while ( my $line = <DBNAMEF> ) {
    chomp $line;
    @temp_array = split( /\s+/, $line );
    $dbname_override{ $temp_array[0] } = $temp_array[1];
  }

  close(DBNAMEF);
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
    if ($dbname_override{$full_db_name}) {
        $full_db_name = $dbname_override{$full_db_name};
    }
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

    
    # check if table exists
    my $test_sql = "select count(1) from information_schema.tables where table_schema = ? and table_name = 'analysis_description'";
    my $test_sth = $dbh->prepare($test_sql);
    $test_sth->execute($dbname);
    my ($table_exists) = $test_sth->fetchrow_array();
    if ($table_exists) {
        my $file_path = $dumppath . "/" . $dbname . 'analysis_description.sql';
        my $file_exists = 0;
        my $response;
        if (-e $file_path) {
                print("file $file_path already exists, overwrite? (y/n)\n");
                $file_exists = 1;
                $response = <>;
                chomp $response;
        }
        if ( !$file_exists or $response eq 'y') {
                open(BKUPFILE, ">$file_path") or die("Failed to open file $file_path for writing\n");   
                my $cmd = "mysqldump -h $host -u $user -p$pass $dbname analysis_description";
                my $result = `$cmd`;
                print BKUPFILE $result;
                close BKUPFILE;
                if ($result !~ /Dump completed/) { 
                        print("back up failed, check file $file_path for details\n");
                        next;
                } else {
                        print("$full_table_name dumped out to file $file_path\n");
                }
        } else {
          print("skipping update for table analysis_description\n");
          next; 
        }
      } else {
          print("table analysis_description does not exist in database $dbname\n");
          next;
      }


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

    }
    #delete the backup table
    $dbh->do(
           sprintf( 'DROP TABLE IF EXISTS %s', $full_table_name_bak ) );

  } continue {
    print("\n");
  }

  $dbh->disconnect();
}

print "To restore a table from dump login to the database and use command: source {dump file name};\n";
