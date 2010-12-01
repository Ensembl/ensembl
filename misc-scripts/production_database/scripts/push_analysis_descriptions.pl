#!/usr/bin/env perl

use strict;
use warnings;

use DBI qw( :sql_types );
use File::Spec::Functions;
use Getopt::Long qw( :config no_ignore_case );
use IO::File;
use POSIX qw( floor ceil );

my $outdir = 'fix-analysis_description';

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
TODO: About info.
ABOUT_END
}

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

my @dbtypes = ( 'core', 'otherfeatures', 'cdna', 'vega' );

my %db_handles;

my %master;

{
  my $dsn = sprintf( 'DBI:mysql:host=%s;port=%d;database=%s',
                     $master, $dbport, 'ensembl_production' );
  my $dbh =
    DBI->connect( $dsn, $dbuser, $dbpass, { 'PrintError' => 1 } );

  my $sth =
    $dbh->prepare(   'SELECT logic_name, description, display_label '
                   . 'FROM analysis_description' );

  $sth->execute();

  my ( $logic_name, $description, $display_label );

  $sth->bind_columns( \( $logic_name, $description, $display_label ) );

  while ( $sth->fetch() ) {
    my $logic_name_lc = lc($logic_name);

    $master{$logic_name_lc} = { 'logic_name'    => $logic_name_lc,
                                'description'   => $description,
                                'display_label' => $display_label };
  }
}

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

      my $sth2 = $dbh->prepare(
                 sprintf(
                   'SELECT a.logic_name, '
                     . 'ad.description, '
                     . 'ad.display_label '
                     . 'FROM %s a JOIN %s ad USING (analysis_id)',
                   $dbh->quote_identifier( undef, $dbname, 'analysis' ),
                   $dbh->quote_identifier(
                                  undef, $dbname, 'analysis_description'
                   ) ) );

      $sth2->execute();

      my ( $logic_name, $description, $display_label );

      $sth2->bind_columns(
                       \( $logic_name, $description, $display_label ) );

      while ( $sth2->fetch() ) {
        my $logic_name_lc = lc($logic_name);

        if ( exists( $master{$logic_name_lc} )
             && $logic_name ne $master{$logic_name_lc}{'logic_name'} )
        {
          # Wrong capitalization in analysis.logic_name.

          display_banner( '-', $dbname );

          printf(
               "==> The logic name '%s' should be all in lower case.\n",
               $logic_name );

          push( @{ $sql{$dbname} },
                sprintf( "-- Updating (lower-casing) logic_name '%s'\n"
                           . "UPDATE %s\n\t"
                           . "SET logic_name = %s\n\t"
                           . "WHERE logic_name = %s;\n\n",
                         $logic_name,
                         $dbh->quote_identifier(
                                              undef, $dbname, 'analysis'
                         ),
                         $dbh->quote( $logic_name_lc, SQL_VARCHAR ),
                         $dbh->quote( $logic_name,    SQL_VARCHAR ) ) );
        }

        if ( !exists( $master{$logic_name_lc} ) ) {
          # Missing in master.

          display_banner( '-', $dbname );

          print("==> Description MISSING IN MASTER:\n");

          printf( "==> logic_name = '%s'\n"
                    . "==> description = '%s'\n"
                    . "==> display_label = '%s'\n",
                  $logic_name, $description, $display_label );

          push( @{ $sql{$dbname} },
                sprintf(
                    "#!! -- MASTER: Inserting logic_name '%s'\n"
                      . "#!! INSERT INTO %s (\n"
                      . "#!! \tlogic_name, description, display_label\n"
                      . "#!! ) VALUES (\n"
                      . "#!! \t%s,\n"
                      . "#!! \t%s,\n"
                      . "#!! \t%s\n"
                      . "#!! );\n\n",
                    $logic_name_lc,
                    $dbh->quote_identifier( undef, 'ensembl_production',
                                            'analysis_description' ),
                    $dbh->quote( $logic_name_lc, SQL_VARCHAR ),
                    $dbh->quote( $description,   SQL_VARCHAR ),
                    $dbh->quote( $display_label, SQL_VARCHAR ) ) );

          $master{$logic_name_lc} = { 'logic_name'    => $logic_name_lc,
                                      'description'   => $description,
                                      'display_label' => $display_label
          };
        } else {
          # Compare all fields.

          if ( (   !defined($description)
                 && defined( $master{$logic_name_lc}{'description'} ) )
               || ( defined($description)
                 && !defined( $master{$logic_name_lc}{'description'} ) )
               || (    defined($description)
                    && defined( $master{$logic_name_lc}{'description'} )
                    && $description ne
                    $master{$logic_name_lc}{'description'} ) )
          {
            # Description differs.
            display_banner( '-', $dbname );

            printf( "==> Description differs for logic_name '%s':\n",
                    $logic_name );
            printf( "==> In table:\t%s\n", $description || 'NULL' );
            printf( "==> In master:\t%s\n",
                    $master{$logic_name_lc}{'description'} );

            push( @{ $sql{$dbname} },
                  sprintf(
                         "-- Updating description for logic_name '%s'\n"
                           . "UPDATE %s ad,\n"
                           . "       %s a\n"
                           . "SET ad.description = %s\n"
                           . "WHERE a.logic_name = %s\n"
                           . "AND ad.analysis_id = a.analysis_id;\n",
                         $logic_name_lc,
                         $dbh->quote_identifier(
                                  undef, $dbname, 'analysis_description'
                         ),
                         $dbh->quote_identifier(
                                              undef, $dbname, 'analysis'
                         ),
                         $dbh->quote(
                                 $master{$logic_name_lc}{'description'},
                                 SQL_VARCHAR
                         ),
                         $dbh->quote( $logic_name_lc, SQL_VARCHAR )
                  ),
                  sprintf( "-- previous value was '%s'\n",
                           $description || 'NULL' ),
                  "\n"
            );

          } ## end if ( ( !defined($description...)))

          if ( (
                 !defined($display_label)
               && defined( $master{$logic_name_lc}{'display_label'} ) )
             || ( defined($display_label)
               && !defined( $master{$logic_name_lc}{'display_label'} ) )
             || (    defined($display_label)
                  && defined( $master{$logic_name_lc}{'display_label'} )
                  && $display_label ne
                  $master{$logic_name_lc}{'display_label'} ) )
          {
            # Display label differs.
            display_banner( '-', $dbname );

            printf( "==> display_label differs for logic_name '%s':\n",
                    $logic_name );
            printf( "==> In table:\t%s\n", $display_label || 'NULL' );
            printf( "==> In master:\t%s\n",
                    $master{$logic_name_lc}{'display_label'} );

            push( @{ $sql{$dbname} },
                  sprintf(
                       "-- Updating display_label for logic_name '%s'\n"
                         . "UPDATE %s ad,\n"
                         . "       %s a\n"
                         . "SET ad.display_label = %s\n"
                         . "WHERE a.logic_name = %s\n"
                         . "AND ad.analysis_id = a.analysis_id;\n",
                       $logic_name_lc,
                       $dbh->quote_identifier(
                                  undef, $dbname, 'analysis_description'
                       ),
                       $dbh->quote_identifier(undef, $dbname, 'analysis'
                       ),
                       $dbh->quote(
                               $master{$logic_name_lc}{'display_label'},
                               SQL_VARCHAR
                       ),
                       $dbh->quote( $logic_name_lc, SQL_VARCHAR )
                  ),
                  sprintf( "-- previous value was '%s'\n",
                           $display_label ),
                  "\n"
            );

          } ## end if ( ( !defined($display_label...)))
        } ## end else [ if ( !exists( $master{...}))]

      } ## end while ( $sth2->fetch() )
    } ## end while ( $sth->fetch() )

  } ## end foreach my $dbtype (@dbtypes)
} ## end foreach my $server (@servers)

if ( scalar( keys(%sql) ) > 0 ) {
  if ( !-d $outdir ) {
    mkdir($outdir);
  }

  foreach my $db_name ( keys(%sql) ) {
    my $filename =
      catfile( $outdir, sprintf( 'fix-%s.sql', $db_name ) );

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
