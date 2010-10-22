#!/usr/bin/env perl

use strict;
use warnings;

use DBI qw( :sql_types );
use Getopt::Long qw( :config no_ignore_case );

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

my @db_handles;

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

  push( @db_handles, $dbh );
}

foreach my $dbh (@db_handles) {
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

          printf( "%s <=> %s\n",
                  $logic_name,
                  $master{ lc($logic_name) }{'logic_name'} );

        } elsif ( !exists( $master{ lc($logic_name) } ) ) {
          # Missing in master.

          printf( "INSERT INTO %s (\n"
                    . "\tlogic_name, description, display_label\n) "
                    . "VALUES (\n\t%s,\n\t%s,\n\t%s\n);\n",
                  $dbh->quote_identifier(
                    undef, sprintf( 'ensembl_production_%d', $release ),
                    'analysis_description' ),
                  $dbh->quote( $logic_name,    SQL_VARCHAR ),
                  $dbh->quote( $description,   SQL_VARCHAR ),
                  $dbh->quote( $display_label, SQL_VARCHAR ) );

          $master{ lc($logic_name) } = {
                                       'logic_name'    => $logic_name,
                                       'description'   => $description,
                                       'display_label' => $display_label
          };
        }

      } ## end while ( $sth2->fetch() )

    } ## end while ( $sth->fetch() )

  } ## end foreach my $dbtype (@dbtypes)
} ## end foreach my $dbh (@db_handles)

END {
  foreach my $dbh (@db_handles) {
    $dbh->disconnect();
  }
}
