#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use DBI qw( :sql_types );
use Getopt::Long qw( :config no_ignore_case );
use IO::File;
use POSIX qw( floor ceil );

sub usage {
  my $padding = ' ' x length($0);

  print <<USAGE_END;
Usage:
  $0 --release NN --master master_server \\
  $padding --server server1 --server server2 [...] \\
  $padding --dbport 3306 --dbuser user --dbpass passwd

or
  $0 --help

or
  $0 --about

where

  --release/-r  The current release (required).

  --master/-m   The master server where the production database lives
                (optional, default is 'ens-staging1').
  --server/-s   A database server (optional, may occur several times,
                default is 'ens-staging1' and 'ens-staging2').

  --dbport/-P   The port to connect to (optional, default is '3306').

  --dbuser/-u   The (read-only) user to connect as (optional,
                default is 'ensro').
  --dbpass/-p   The password to connect with (optional, no default).

  --help/-h     Displays this help text.

  --about/-a    Display a text about this program (what it does etc.).
USAGE_END
} ## end sub usage

sub about {
  print <<ABOUT_END;
Run the program with --help to get information about available command
line switches.

This program takes the master tables from the production database
and compares it to the corresponding tables on the given servers (by
default, the staging servers).

The program will display any discrepancies on the display
while writing SQL to files in the current directory that will
correct the discrepancies.

Each SQL patch file will have the generic name "fix-DBNAME.sql"
where "DBNAME" is the name of the database, e.g.,
"fix-oryctolagus_cuniculus_otherfeatures_60_3.sql".

A discrepancy is patched by

1)  Insertion into the master table in the production database in the
    case where a new entry has been added to a database without being
    added to the master table.

2)  Insertion into the database table in the case where a new master
    entry is missing in the database.

3)  Updating the database entry in the case where an entry (identified by
    its primary key only) differs in any of its fields.

The SQL patch files may then be used to patch the databases:

  \$ mysql -h server -u user -ppass < fix-DBNAME.sql




                    BE SURE TO REVIEW THESE SQL PATCH FILES
                        (along with the program output)

                            WITH YOUR EYE AND BRAIN

                              BEFORE APPLYING THEM


ABOUT_END
} ## end sub about

sub fetch_table {
  my ( $dbh, $dbname, $table ) = @_;

  my $sth = $dbh->prepare(
             sprintf( 'SELECT * FROM %s',
                      $dbh->quote_identifier( undef, $dbname, $table ) )
  );

  $sth->execute();

  my %table_hash;

  $table =~ s/^master_//;

  my $pk = sprintf( '%s_id', $table );
  while ( my $row = $sth->fetchrow_hashref() ) {
    if ( !exists( $row->{$pk} ) ) {
      die( sprintf( "Can not find expected primary key '%s'", $pk ) );
    }

    $table_hash{ $row->{$pk} } = $row;
  }

  return \%table_hash;
} ## end sub fetch_table

sub display_banner {
  my ( $char, $text ) = @_;

  printf( "%s %s %s\n",
          $char x ( 39 - floor( length($text)/2 ) ),
          $text, $char x ( 39 - ceil( length($text)/2 ) ) );
}

my $release;
my @servers = ( 'ens-staging1', 'ens-staging2' );
my $master = 'ens-staging1';

my $dbport = '3306';
my ( $dbuser, $dbpass ) = ( 'ensro', undef );

my $opt_help  = 0;
my $opt_about = 0;

if ( !GetOptions( 'release|r=i' => \$release,
                  'master|m=s'  => \$master,
                  'server|s=s@' => \@servers,
                  'dbuser|u=s'  => \$dbuser,
                  'dbpass|p=s'  => \$dbpass,
                  'dbport|P=s'  => \$dbport,
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

my @tables =
  ( 'attrib_type', 'external_db', 'misc_set', 'unmapped_reason' );
my @dbtypes = ( 'core', 'otherfeatures', 'cdna', 'vega' );

my %master;
{
  my $dsn = sprintf( 'DBI:mysql:host=%s;port=%d', $master, $dbport );
  my $dbh =
    DBI->connect( $dsn, $dbuser, $dbpass, { 'PrintError' => 1 } );

  foreach my $table (@tables) {
    my $master_table = sprintf( 'master_%s', $table );
    $master{$table} = fetch_table( $dbh,
                                   sprintf( 'ensembl_production_%d',
                                            $release ),
                                   $master_table );
  }
}

my %db_handles;
foreach my $server (@servers) {
  my $dsn = sprintf( 'DBI:mysql:host=%s;port=%d', $server, $dbport );
  my $dbh =
    DBI->connect( $dsn, $dbuser, $dbpass, { 'PrintError' => 1 } );

  $db_handles{$server} = $dbh;
}

my %sql;
foreach my $server (@servers) {
  printf( "###> Looking at '%s'\n", $server );
  my $dbh = $db_handles{$server};

  my $sth = $dbh->prepare('SHOW DATABASES LIKE ?');

  foreach my $dbtype (@dbtypes) {
    $sth->bind_param( 1,
                      sprintf( '%%\\_%s\\_%d\\_%%', $dbtype, $release ),
                      SQL_VARCHAR );

    $sth->execute();

    my $dbname;
    $sth->bind_col( 1, \$dbname );

    while ( $sth->fetch() ) {
      printf( "##> Processing '%s'\n", $dbname );

      foreach my $table (@tables) {
        my $csth = $dbh->column_info( undef, $dbname, $table, '%' );
        my $colinfo = $csth->fetchall_hashref( ['COLUMN_NAME'] );

        my %table = %{ fetch_table( $dbh, $dbname, $table ) };

        foreach
          my $pk ( sort { $a <=> $b } keys( %{ $master{$table} } ) )
        {
          if ( !exists( $table{$pk} ) ) {
            my $row    = $master{$table}{$pk};
            my @fields = sort( keys( %{$row} ) );

            push(
              @{ $sql{$dbname} },
              sprintf( "-- insert %s_id=%d in %s\n",
                       $table, $pk, $table ),
              sprintf(
                "INSERT INTO %s (\n\t%s\n) VALUES (\n\t%s\n);\n",
                $dbh->quote_identifier( undef, $dbname, $table ),
                join( ",\n\t",
                      map { $dbh->quote_identifier($_) } @fields ),
                join(
                  ",\n\t",
                  map {
                    $dbh->quote( $row->{$_},
                                 $colinfo->{$_}{'DATA_TYPE'} )
                    } @fields ) ) );

          }
        } ## end foreach my $pk ( sort { $a ...})

        foreach my $pk ( sort { $a <=> $b } keys(%table) ) {
          my $master_row = $master{$table}{$pk};
          my $row        = $table{$pk};

          if ( $pk == 0 ) {
            display_banner( '-', sprintf( '%s.%s', $dbname, $table ) );

            print( "==> Primary key is ZERO "
                     . "for the following row in DATABASE:\n",
                   Dumper($row),
                   "\n" );
          } else {
            my @fields = sort( keys( %{$row} ) );

            if ( !defined($master_row) ) {
              display_banner( '=',
                              sprintf( '%s.%s', $dbname, $table ) );

              # Find other row in master table that is the same as
              # database table row, but with different primary key.

              my $is_missing = 1;

              foreach my $master_pk ( keys( %{ $master{$table} } ) ) {
                my $master_row = $master{$table}{$master_pk};
                my $is_same    = 1;

                foreach my $field ( sort( keys( %{$master_row} ) ) ) {
                  if ( $field eq sprintf( '%s_id', $table ) ) {
                    # Skip the primary key.
                    next;
                  }

                  if ( $master_row->{$field} ne $row->{$field} ) {
                    $is_same = 0;
                    last;
                  }
                }

                if ($is_same) {
                  printf( "==> Entry with primary key %d "
                          . "is same as entry with primary key %d:\n%s",
                        $pk, $master_pk, Dumper($master_row) );

                  push( @{ $sql{$dbname} },
                        sprintf(
                               "-- Entries with %s_id = %d "
                                 . "should change this to %d\n"
                                 . "-- Useful SQL:\n"
                                 . "-- UPDATE <table> "
                                 . "SET %s_id = %d WHERE %s_id = %s;\n",
                               $table,     $pk,    $master_pk, $table,
                               $master_pk, $table, $pk ) );

                  $is_missing = 0;
                }
              } ## end foreach my $master_pk ( keys...)

              if ($is_missing) {
                print( "==> The following row is MISSING IN MASTER:\n",
                       Dumper($row) );

                push(
                  @{ $sql{$dbname} },
                  sprintf( "#HEADS_UP!# -- MASTER: insert from %s.%s\n",
                           $dbname, $table ),
                  sprintf(
                    "#HEADS_UP!# INSERT INTO %s (\n\t%s\n) "
                      . "VALUES (\n\t%s\n);\n",
                    $dbh->quote_identifier(
                           undef,
                           sprintf( 'ensembl_production_%d', $release ),
                           sprintf( 'master_%s',             $table ) ),
                    join( ",\n#HEADS_UP!# \t",
                          map { $dbh->quote_identifier($_) } @fields ),
                    join(
                      ",\n#HEADS_UP!# \t",
                      map {
                        $dbh->quote( $row->{$_},
                                     $colinfo->{$_}{'DATA_TYPE'} )
                        } @fields ) ) );

                print("\n");
              } ## end if ($is_missing)
            } else {
              my %diff_fields;

              foreach my $field (@fields) {
                if (    defined( $master_row->{$field} )
                     || defined( $row->{$field} ) )
                {
                  if ( (   !defined( $master_row->{$field} )
                         && defined( $row->{$field} ) )
                       || ( defined( $master_row->{$field} )
                            && !defined( $row->{$field} ) )
                       || ( $master_row->{$field} ne $row->{$field} ) )
                  {
                    if ( !(    $table eq 'external_db'
                            && $field eq 'db_release' ) )
                    {
                      $diff_fields{$field} = $master_row->{$field};
                    }
                  }
                }
              }

              if ( scalar( keys(%diff_fields) ) > 0 ) {
                display_banner( '=',
                                sprintf( '%s.%s', $dbname, $table ) );

                # Find other row in master table that is the same as
                # database table row, but with different primary key.

                my $is_missing = 1;

                foreach my $master_pk ( keys( %{ $master{$table} } ) ) {
                  my $master_row = $master{$table}{$master_pk};
                  my $is_same    = 1;

                  foreach my $field ( sort( keys( %{$master_row} ) ) ) {
                    if ( $field eq sprintf( '%s_id', $table ) ) {
                      # Skip the primary key.
                      next;
                    }

                    if ( $master_row->{$field} ne $row->{$field} ) {
                      $is_same = 0;
                      last;
                    }
                  }

                  if ($is_same) {
                    printf( "==> Entry with primary key %d "
                          . "is same as entry with primary key %d:\n%s",
                        $pk, $master_pk, Dumper($master_row) );

                    push( @{ $sql{$dbname} },
                          sprintf(
                               "-- Entries with %s_id = %d "
                                 . "should change this to %d\n"
                                 . "-- Useful SQL:\n"
                                 . "-- UPDATE <table> "
                                 . "SET %s_id = %d WHERE %s_id = %s;\n",
                               $table,     $pk,    $master_pk, $table,
                               $master_pk, $table, $pk ) );

                    $is_missing = 0;
                  }
                } ## end foreach my $master_pk ( keys...)

                if ($is_missing) {
                  printf( "==> The following row differs in %s.\n",
                          join( ', ', keys(%diff_fields) ) );
                  print( "==> MASTER row:\n",   Dumper($master_row),
                         "==> DATABASE row:\n", Dumper($row) );

                  push(
                    @{ $sql{$dbname} },
                    sprintf( "-- update %s in %s\n",
                             join( ', ', keys(%diff_fields) ), $table ),
                    sprintf(
                      "UPDATE %s\nSET %s\nWHERE %s_id = %d;\n",
                      $dbh->quote_identifier( undef, $dbname, $table ),
                      join(
                        ', ',
                        map {
                          sprintf( '%s = %s',
                                   $_,
                                   $dbh->quote( $diff_fields{$_} ),
                                   $colinfo->{$_}{'DATA_TYPE'} )
                          }
                          keys(%diff_fields) ),
                      $table,
                      $pk ) );

                  print("\n");
                } ## end if ($is_missing)

              } ## end if ( scalar( keys(%diff_fields...)))
            } ## end else [ if ( !defined($master_row...))]

          } ## end else [ if ( $pk == 0 ) ]

        } ## end foreach my $pk ( sort { $a ...})
      } ## end foreach my $table (@tables)

    } ## end while ( $sth->fetch() )
  } ## end foreach my $dbtype (@dbtypes)

} ## end foreach my $dbh (@db_handles)

if ( scalar( keys(%sql) ) > 0 ) {
  foreach my $db_name ( keys(%sql) ) {
    my $filename = sprintf( 'fix-%s.sql', $db_name );
    printf( "==> Writing SQL to '%s'\n", $filename );
    my $out = IO::File->new( $filename, 'w' );
    $out->print( @{ $sql{$db_name} } );
    $out->close();
  }
} else {
  print("Nothing to do, all seems ok\n");
}

END {
  foreach my $dbh ( values(%db_handles) ) {
    $dbh->disconnect();
  }
}
