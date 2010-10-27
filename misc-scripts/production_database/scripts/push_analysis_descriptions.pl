#!/usr/bin/env perl

use strict;
use warnings;

use DBI qw( :sql_types );
use File::Spec::Functions;
use Getopt::Long qw( :config no_ignore_case );
use IO::File;
use POSIX qw( floor ceil );

sub usage {
  print <<USAGE_END;
TODO: Usage info.
USAGE_END
}

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
                     $master, $dbport,
                     sprintf( 'ensembl_production_%d', $release ) );
  my $dbh =
    DBI->connect( $dsn, $dbuser, $dbpass, { 'PrintError' => 1 } );

  my $sth =
    $dbh->prepare(   'SELECT logic_name, description, display_label '
                   . 'FROM analysis_description' );

  $sth->execute();

  my ( $logic_name, $description, $display_label );

  $sth->bind_columns( \( $logic_name, $description, $display_label ) );

  while ( $sth->fetch() ) {
    $master{ lc($logic_name) } = { 'logic_name'    => $logic_name,
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
        if ( exists( $master{ lc($logic_name) } )
            && $logic_name ne $master{ lc($logic_name) }{'logic_name'} )
        {
          # Wrong capitalization in analysis.logic_name.

          display_banner( '-', $dbname );

          printf(
               "==> The logic name '%s' should be all in lower case.\n",
               $logic_name );

          push( @{ $sql{$dbname} },
                sprintf( "-- Updating logic_name '%s'\n"
                           . "UPDATE %s\n\t"
                           . "SET logic_name = %s\n\t"
                           . "WHERE logic_name = %s;\n\n",
                         lc($logic_name),
                         $dbh->quote_identifier(
                                              undef, $dbname, 'analysis'
                         ),
                         $dbh->quote( lc($logic_name), SQL_VARCHAR ),
                         $dbh->quote( $logic_name,     SQL_VARCHAR ) )
          );
        }

        if ( !exists( $master{ lc($logic_name) } ) ) {
          # Missing in master.

          display_banner( '-', $dbname );

          print("==> Description MISSING IN MASTER:\n");

          printf( "==> logic_name = '%s'\n"
                    . "==> description = '%s'\n"
                    . "==> display_label = '%s'\n",
                  $logic_name, $description, $display_label );

          push( @{ $sql{$dbname} },
                sprintf(
                      "-- Inserting logic_name '%s' in master\n"
                        . "INSERT INTO %s (\n"
                        . "\tlogic_name, description, display_label\n) "
                        . "VALUES (\n\t%s,\n\t%s,\n\t%s\n);\n\n",
                      lc($logic_name),
                      $dbh->quote_identifier(
                           undef,
                           sprintf( 'ensembl_production_%d', $release ),
                           'analysis_description' ),
                      $dbh->quote( lc($logic_name), SQL_VARCHAR ),
                      $dbh->quote( $description,    SQL_VARCHAR ),
                      $dbh->quote( $display_label,  SQL_VARCHAR ) ) );

          $master{ lc($logic_name) } = {
                                       'logic_name'  => lc($logic_name),
                                       'description' => $description,
                                       'display_label' => $display_label
          };
        } else {
          # Compare all fields.

          if (
             $description ne $master{ lc($logic_name) }{'description'} )
          {
            # Description differs.
            display_banner( '-', $dbname );

            printf( "==> Description differs for logic_name '%s':\n",
                    $logic_name );
            printf( "==> In table:\t%s\n", $description );
            printf( "==> In master:\t%s\n",
                    $master{ lc($logic_name) }{'description'} );

            push( @{ $sql{$dbname} },
                  sprintf(
                        "-- Updating description for logic_name '%s'\n"
                          . "UPDATE %s ad, %s a\n\t"
                          . "SET ad.description = %s\n\t"
                          . "WHERE a.logic_name = %s\n\t"
                          . "AND ad.analysis_id = a.analysis_id;\n\n",
                        lc($logic_name),
                        $dbh->quote_identifier(
                                  undef, $dbname, 'analysis_description'
                        ),
                        $dbh->quote_identifier(
                                              undef, $dbname, 'analysis'
                        ),
                        $dbh->quote(
                              $master{ lc($logic_name) }{'description'},
                              SQL_VARCHAR ),
                        $dbh->quote( lc($logic_name), SQL_VARCHAR ) ) );
          } ## end if ( $description ne $master...)

          if ( $display_label ne
               $master{ lc($logic_name) }{'display_label'} )
          {
            # Display label differs.
            display_banner( '-', $dbname );

            printf( "==> display_label differs for logic_name '%s':\n",
                    $logic_name );
            printf( "==> In table:\t%s\n", $display_label );
            printf( "==> In master:\t%s\n",
                    $master{ lc($logic_name) }{'display_label'} );

            push( @{ $sql{$dbname} },
                  sprintf(
                      "-- Updating display_label for logic_name '%s'\n"
                        . "UPDATE %s ad, %s a\n\t"
                        . "SET ad.display_label = %s\n\t"
                        . "WHERE a.logic_name = %s\n\t"
                        . "AND ad.analysis_id = a.analysis_id;\n\n",
                      lc($logic_name),
                      $dbh->quote_identifier(
                                  undef, $dbname, 'analysis_description'
                      ),
                      $dbh->quote_identifier( undef, $dbname, 'analysis'
                      ),
                      $dbh->quote(
                            $master{ lc($logic_name) }{'display_label'},
                            SQL_VARCHAR ),
                      $dbh->quote( lc($logic_name), SQL_VARCHAR ) ) );
          } ## end if ( $display_label ne...)
        } ## end else [ if ( !exists( $master{...}))]

      } ## end while ( $sth2->fetch() )
    } ## end while ( $sth->fetch() )

  } ## end foreach my $dbtype (@dbtypes)
} ## end foreach my $server (@servers)

my $outdir = 'analysis_desc_fixes';
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
