#!/usr/bin/env perl

use strict;
use warnings;

use DBI qw( :sql_types );
use Getopt::Long qw( :config no_ignore_case );

local $| = 1;

my %master_tables = ( 'attrib_type'     => 1,
                      'external_db'     => 1,
                      'misc_set'        => 1,
                      'unmapped_reason' => 1 );
my @tables;

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
my $cleanup  = 0;
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
                  'table|t=s'      => \@tables,
                  'verbose|v!'     => \$verbose,
                  'core=i'         => \$core,
		  'dumppath|dp=s'   => \$dumppath)
     || !(
           defined($host)
        && defined($user)
        && defined($pass)
        && ( defined($dbname) || defined($dbpattern) || defined($core) )
        && defined($mhost)
        && defined($muser) 
        && defined($dumppath)) )
{
  my $indent = ' ' x length($0);
  print <<USAGE_END;
This script copies tables 'attrib_type', 'external_db', 'misc_set'
and 'unmapped_reason' from the production database into a user-defined
database.

Usage:

  $0 -h host [-P port] \\
  $indent -u user [-p password]
  $indent -d database | --pattern pattern \\
  $indent -dp dumppath
  $indent [-mh host] [-mP port] \\
  $indent [-mu user] [-mp password] [-md database] \\
  $indent [-t table] [-t table] [-t ...] \\
  $indent [-v]

  -h / --host       User database server host
  -P / --port       User database server port (optional, default is 3306)

  -u / --user       User username (must have write-access)
  -p / --pass       User password

  -d / --database   User database name or SQL pattern
                    e.g. --database="homo_sapiens_rnaseq_62_37g"
                    or   --database="%core_62%"

  -dp / --dumppath  Dumpout path. Dump out tables into the specified directory path

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

  -t / --table      A specific table to update, may occur several times,
                    must be one of the tables attrib_type, external_db,
                    misc_set, or unmapped_reason

  -v / --verbose    Be verbose, display every SQL statement as they
                    are executed (on standard error)


USAGE_END

  die(   "Need the following options: "
       . "-h -u -p and -d (or --pattern)\n" );

} ## end if ( !GetOptions( 'mhost|mh=s'...))

if (@tables) {
  foreach my $table (@tables) {
    if ( !exists( $master_tables{$table} ) ) {
      die( sprintf( "Invalid table specified: '%s'\n", $table ) );
    }
  }
} else {
  @tables = keys(%master_tables);
}

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

  foreach my $table (@tables) {
    my $sth = $dbh->prepare(
                   sprintf( 'SELECT * FROM %s',
                     $dbh->quote_identifier( undef, $mdbname, $table ) )
    );

    $sth->execute();

    while ( my $row = $sth->fetchrow_arrayref() ) {
      push( @{ $data{$table} }, [ @{$row} ] );
    }
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

    print( '=' x 80, "\n" );
    printf( "\t%s\n", $dbname );
    print( '=' x 80, "\n" );

    foreach my $table ( keys(%data) ) {
      printf( "==> %s: ", $table );

      my $full_table_name =
        $dbh->quote_identifier( undef, $dbname, $table );
      my $full_table_name_bak =
        $dbh->quote_identifier( undef, $dbname, $table . '_bak' );
      my $key_name = $table . '_id';


      # check if table exists
      my $test_sql = "select count(1) from information_schema.tables where table_schema = ? and table_name = ?";
      my $test_sth = $dbh->prepare($test_sql);
      $test_sth->execute($dbname, $table);
      my ($table_exists) = $test_sth->fetchrow_array();
      if ($table_exists) {
      	my $file_path = $dumppath . "/" . $dbname . $table .'.sql';
      	open(BKUPFILE, ">$file_path") or die("Failed to open file $file_path for writing\n");   
      	my $cmd = "mysqldump -h $host -u $user -p$pass $dbname $table";
      	my $result = `$cmd`;
      	print BKUPFILE $result;
      	close BKUPFILE;
      	if ($result !~ /Dump completed/) { 
	  print("back up failed, check file $file_path for details\n");
	  next;
      	} else {
	  print("$full_table_name_bak dumped out to file $file_path\n");
      	}
      } else {
	  print("table $table does not exist in database $dbname\n");
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

      # Get column information.
      my $colinfo_sth =
        $dbh->column_info( undef, $dbname, $table, '%' );
      my $colinfo =
        $colinfo_sth->fetchall_hashref( ['ORDINAL_POSITION'] );

      my $numcols = scalar( keys( %{$colinfo} ) );

      # For each row read from the master table,
      # issue an INSERT statement.
      foreach my $row ( @{ $data{$table} } ) {
        my $insert_statement = sprintf(
          'INSERT INTO %s (%s) VALUES (%s)',
          $full_table_name,
          join( ', ',
                map { $colinfo->{$_}{'COLUMN_NAME'} } 1 .. $numcols ),
          join(
            ', ',
            map {
              $dbh->quote( $row->[ $_ - 1 ],
                           $colinfo->{$_}{'DATA_TYPE'} )
              } ( 1 .. $numcols ) ) );

        if ($verbose) {
          printf( STDERR "EXECUTING: %s\n", $insert_statement );
        }
        $dbh->do($insert_statement);
      }

      print("<inserted data>");

      {
        my $statement = sprintf( 'SELECT %s '
                                   . 'FROM %s '
                                   . 'LEFT JOIN %s t USING (%s) '
                                   . 'WHERE t.%s IS NULL '
                                   . 'ORDER BY %s',
                                 $key_name,
                                 $full_table_name,
                                 $full_table_name_bak,
                                 $key_name,
                                 $key_name,
                                 $key_name );

        my $sth2 = $dbh->prepare($statement);

        if ($verbose) {
          printf( STDERR "EXECUTING: %s\n", $statement );
        }
        $sth2->execute();

        my $key;
        $sth2->bind_col( 1, \$key );

        my @keys;
        while ( $sth2->fetch() ) {
          push( @keys, $key );
        }

        if (@keys) {
          print("\n");
          print("New data inserted:\n");
          printf( "SELECT * FROM %s WHERE %s_id IN (%s);\n",
                  $table, $table, join( ',', @keys ) );
          print("\n");
        }
      }
      {
        my $statement = sprintf( 'SELECT %s '
                                   . 'FROM %s '
                                   . 'LEFT JOIN %s t USING (%s) '
                                   . 'WHERE t.%s IS NULL '
                                   . 'ORDER BY %s',
                                 $key_name,
                                 $full_table_name_bak,
                                 $full_table_name,
                                 $key_name,
                                 $key_name,
                                 $key_name );

        my $sth2 = $dbh->prepare($statement);

        if ($verbose) {
          printf( STDERR "EXECUTING: %s\n", $statement );
        }
        $sth2->execute();

        my $key;
        $sth2->bind_col( 1, \$key );

        my @keys;
        while ( $sth2->fetch() ) {
          push( @keys, $key );
        }

        if (@keys) {
          print("\n");
          print( '-' x 40, "\n" );
          print("!! Old data deleted:\n");
          printf( "SELECT * FROM %s WHERE %s_id IN (%s);\n",
                  $table . '_bak',
                  $table, join( ',', @keys ) );
          print( '-' x 40, "\n" );
          print("\n");
        }
      }
      # delete the backup table  
      $dbh->do(
           sprintf( 'DROP TABLE IF EXISTS %s', $full_table_name_bak ) );
  

    } continue {
      print("\n");
    }

    print("\n");

  } ## end while ( $sth->fetch() )

  $dbh->disconnect();
}
