#!/usr/bin/env perl

use strict;
use warnings;

use Data::Dumper;
use DBI qw( :sql_types );
use File::Spec::Functions;
use Getopt::Long qw( :config no_ignore_case );
use IO::File;
use POSIX qw( floor ceil );

$Data::Dumper::Terse    = 1;
$Data::Dumper::Useqq    = 0;
$Data::Dumper::Indent   = 0;
$Data::Dumper::Deparse  = 0;
$Data::Dumper::Sortkeys = 1;

my $outdir = 'fix-analysis_description';

sub usage {
  my $padding = ' ' x length($0);

  print <<USAGE_END;

Usage:

  $0 --release NN --master master_server[:port] \\
  $padding --server server1[:port1] --server server2[:port2] [...] \\
  $padding --dbuser user --dbpass passwd

or
  $0 --help

or
  $0 --about

where

  --release/-r  The current release (required).

  --master/-m   The master server where the production database lives
                (optional, default is 'ens-staging1').  Specifying the
                port is optional, the default port is 3306.

  --server/-s   A database server (optional, may occur several times,
                default is 'ens-staging1' and 'ens-staging2').
                Specifying the port is optional, the default port is
                3306.

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
my @servers;
my $master = 'ens-staging1';

my ( $dbuser, $dbpass ) = ( 'ensro', undef );

my $opt_help  = 0;
my $opt_about = 0;

if ( !GetOptions( 'release|r=i' => \$release,
                  'master|m=s'  => \$master,
                  'server|s=s@' => \@servers,
                  'dbuser|u=s'  => \$dbuser,
                  'dbpass|p=s'  => \$dbpass,
                  'help|h!'     => \$opt_help,
                  'about!'      => \$opt_about
     )
     || $opt_help
  )
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

if ( !@servers ) {
  @servers = ( 'ens-staging1', 'ens-staging2' );
}

my @dbtypes = ( 'core', 'otherfeatures', 'cdna', 'vega', 'rnaseq' );

my %db_handles;

my %master;
my %master_r;

{
  my ( $host, $port ) = ( $master =~ /^([^:]+):?(\d+)?/ );

  if ( !defined($port) ) {
    $port = 3306;
  }

  my $dsn = sprintf( 'DBI:mysql:host=%s;port=%d;database=%s',
                     $host, $port, 'ensembl_production' );
  my $dbh =
    DBI->connect( $dsn, $dbuser, $dbpass, { 'PrintError' => 1 } );

  my $sth = $dbh->prepare(
    q(
SELECT  full_db_name,
        logic_name,
        description,
        display_label,
        displayable,
        web_data
FROM full_analysis_description
) );

  $sth->execute();

  my ( $full_db_name,  $logic_name,  $description,
       $display_label, $displayable, $web_data );

  $sth->bind_columns( \( $full_db_name,  $logic_name,  $description,
                         $display_label, $displayable, $web_data
                      ) );

  while ( $sth->fetch() ) {
    my $logic_name_lc = lc($logic_name);

    $master{$full_db_name}{$logic_name_lc} = {
                                      'full_db_name'  => $full_db_name,
                                      'logic_name'    => $logic_name_lc,
                                      'description'   => $description,
                                      'display_label' => $display_label,
                                      'displayable'   => $displayable,
                                      'web_data'      => $web_data };

    $master_r{$logic_name_lc}{$full_db_name} =
      $master{$full_db_name}{$logic_name_lc};
  }

}

foreach my $server (@servers) {
  my ( $host, $port ) = ( $server =~ /^([^:]+):?(\d+)?/ );

  if ( !defined($port) ) {
    $port = 3306;
  }

  my $dsn = sprintf( 'DBI:mysql:host=%s;port=%d', $host, $port );
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
      if ( !exists( $master{$dbname} ) ) {
        printf( "!!> Unknown database '%s'\n", $dbname );
        next;
      }

      my $master = $master{$dbname};

      printf( "##> Processing '%s'\n", $dbname );

      my $sth2 = $dbh->prepare(
                 sprintf(
                   'SELECT ad.analysis_id, '
                     . 'a.analysis_id, '
                     . 'a.logic_name,'
                     . 'ad.description,'
                     . 'ad.display_label,'
                     . 'ad.displayable,'
                     . 'ad.web_data '
                     . 'FROM %s a LEFT JOIN %s ad USING (analysis_id)',
                   $dbh->quote_identifier( undef, $dbname, 'analysis' ),
                   $dbh->quote_identifier(
                                  undef, $dbname, 'analysis_description'
                   ) ) );

      $sth2->execute();

      my %row = ( 'exists'        => undef,
                  'analysis_id'   => undef,
                  'logic_name'    => undef,
                  'description'   => undef,
                  'display_label' => undef,
                  'displayable'   => undef,
                  'web_data'      => undef
      );

      $sth2->bind_columns(\( $row{'exists'},        $row{'analysis_id'},
                             $row{'logic_name'},    $row{'description'},
                             $row{'display_label'}, $row{'displayable'},
                             $row{'web_data'} ) );

      while ( $sth2->fetch() ) {
        my $logic_name_lc = lc( $row{'logic_name'} );

        if ( !defined( $row{'exists'} ) ) {
          # An analysis exists that does not have a corresponding
          # analysis_description.

          display_banner( '-', $dbname );

          printf( "==> Analysis '%s' "
                    . "is missing analysis_description entry.\n",
                  $row{'logic_name'} );

          if ( !exists( $master->{$logic_name_lc} ) ) {
            if ( exists( $master_r{$logic_name_lc} ) ) {
              push( @{ $sql{$dbname} },
                    sprintf( "-- WARNING: missing analysis_desciption "
                               . "for logic_name '%s'\n\n",
                             $logic_name_lc ) );
            } else {
              push( @{ $sql{$dbname} },
                    sprintf( "-- WARNING: unknown analysis_desciption "
                               . "for logic_name '%s'\n\n",
                             $logic_name_lc ) );
            }
          } else {

            push(
              @{ $sql{$dbname} },
              sprintf(
                "-- WARNING: Inserting minimal missing "
                  . "analysis_description for '%s'\n"
                  . "-- web_data and displayable fields may be wrong!\n"
                  . "INSERT INTO %s (\n"
                  . "\tanalysis_id,\n"
                  . "\tdescription,\n"
                  . "\tdisplay_label,\n"
                  . "\tdisplayable,\n"
                  . "\tweb_data\n"
                  . ") VALUES (\n"
                  . "\t%s,\n"
                  . "\t%s,\n"
                  . "\t%s,\n"
                  . "\t%s,\n"
                  . "\t%s\n"
                  . ");\n\n",
                $logic_name_lc,
                $dbh->quote_identifier(
                                  undef, $dbname, 'analysis_description'
                ),
                $dbh->quote( $row{'analysis_id'}, SQL_INTEGER ),
                $dbh->quote( $master->{$logic_name_lc}{'description'},
                             SQL_VARCHAR ),
                $dbh->quote( $master->{$logic_name_lc}{'display_label'},
                             SQL_VARCHAR ),
                $dbh->quote( $master->{$logic_name_lc}{'displayable'},
                             SQL_INTEGER ),
                $dbh->quote( $master->{$logic_name_lc}{'web_data'},
                             SQL_VARCHAR ) ) );

          } ## end else [ if ( !exists( $master->...))]

          next;
        } ## end if ( !defined( $row{'exists'...}))

        if ( exists( $master->{$logic_name_lc} )
             && $row{'logic_name'} ne
             $master->{$logic_name_lc}{'logic_name'} )
        {
          # Wrong capitalization in analysis.logic_name.

          display_banner( '-', $dbname );

          printf(
               "==> The logic name '%s' should be all in lower case.\n",
               $row{'logic_name'} );

          push( @{ $sql{$dbname} },
                sprintf( "-- Updating (lower-casing) logic_name '%s'\n"
                           . "UPDATE %s\n"
                           . "SET logic_name = %s\n"
                           . "WHERE logic_name = %s;\n\n",
                         $row{'logic_name'},
                         $dbh->quote_identifier(
                                              undef, $dbname, 'analysis'
                         ),
                         $dbh->quote( $logic_name_lc,     SQL_VARCHAR ),
                         $dbh->quote( $row{'logic_name'}, SQL_VARCHAR )
                ) );
        } ## end if ( exists( $master->...))

        if ( !exists( $master->{$logic_name_lc} ) ) {
          # Missing in master.

          display_banner( '-', $dbname );

          print("==> Analysis description MISSING IN MASTER:\n");

          printf( "==> logic_name = '%s'\n"
                    . "==> description = '%s'\n"
                    . "==> display_label = '%s'\n"
                    . "==> displayable = '%s'\n"
                    . "==> web_data = '%s'\n",
                  $row{'logic_name'}    || 'NULL',
                  $row{'description'}   || 'NULL',
                  $row{'display_label'} || 'NULL',
                  $row{'displayable'}   || 'NULL',
                  $row{'web_data'}      || 'NULL' );

          push( @{ $sql{$dbname} },
                sprintf(
                    "#!! WARNING: THE FOLLOWING IS ONLY A GUIDELINE\n"),
                sprintf(
                    "#!! -- MASTER: Inserting logic_name '%s' "
                      . "into analysis_description\n"
                      . "#!! INSERT INTO %s (\n"
                      . "#!! \tlogic_name, description, display_label\n"
                      . "#!! ) VALUES (\n"
                      . "#!! \t%s,\n"
                      . "#!! \t%s,\n"
                      . "#!! \t%s\n"
                      . "#!! );\n",
                    $logic_name_lc,
                    $dbh->quote_identifier( undef, 'ensembl_production',
                                            'analysis_description' ),
                    $dbh->quote( $logic_name_lc,        SQL_VARCHAR ),
                    $dbh->quote( $row{'description'},   SQL_VARCHAR ),
                    $dbh->quote( $row{'display_label'}, SQL_VARCHAR ) ),
                sprintf( "#!! -- MASTER: Inserting above logic_name "
                           . "into web_data\n"
                           . "#!! INSERT INTO %s (\n"
                           . "#!!\tdata\n"
                           . "#!! ) VALUES (\n"
                           . "#!!\t%s\n"
                           . "#!! );\n",
                         $dbh->quote_identifier(
                                 undef, 'ensembl_production', 'web_data'
                         ),
                         $dbh->quote( (
                                    defined( $row{'web_data'} )
                                    ? Dumper( eval( $row{'web_data'} ) )
                                    : undef ),
                                  SQL_VARCHAR ) ),
                sprintf(
                  "#!! NOW CONNECT THESE IN analysis_web_data! :-)\n\n")
          );

          $master->{$logic_name_lc} = {
                               'logic_name'    => $logic_name_lc,
                               'description'   => $row{'description'},
                               'display_label' => $row{'display_label'},
                               'displayable'   => $row{'displayable'},
                               'web_data'      => $row{'web_data'} };

        } else {
          # Compare all fields.

          my @differs;

          foreach my $field (
                     qw(description display_label displayable web_data))
          {

            if ( $field eq 'web_data' ) {
              if ( (
                     !defined( $row{$field} )
                   && defined( $master->{$logic_name_lc}{$field} ) )
                 || ( defined( $row{$field} )
                      && !defined( $master->{$logic_name_lc}{$field} ) )
                 || ( defined( $row{$field} )
                   && defined( $master->{$logic_name_lc}{$field} )
                   && Dumper( eval( $row{$field} ) ) ne
                   Dumper( eval( $master->{$logic_name_lc}{$field} ) ) )
                )
              {
                # Some field differs.
                push( @differs, $field );
              }
            } else {
              if ( (
                    !defined( $row{$field} )
                  && defined( $master->{$logic_name_lc}{$field} ) )
                || ( defined( $row{$field} )
                     && !defined( $master->{$logic_name_lc}{$field} ) )
                || ( defined( $row{$field} )
                  && defined( $master->{$logic_name_lc}{$field} )
                  && $row{$field} ne $master->{$logic_name_lc}{$field} )
                )
              {
                # Some field differs.
                push( @differs, $field );
              }
            }
          } ## end foreach my $field (...)

          if ( scalar(@differs) > 0 ) {
            my $csth =
              $dbh->column_info( undef, $dbname, 'analysis_description',
                                 '%' );
            my $colinfo = $csth->fetchall_hashref( ['COLUMN_NAME'] );

            display_banner( '-', $dbname );

            $Data::Dumper::Terse  = 0;
            $Data::Dumper::Indent = 1;
            $Data::Dumper::Useqq  = 1;
            printf( "==> Analysis description differs "
                      . "for logic_name '%s'\n",
                    $row{'logic_name'} );
            printf( "==> in the following fields: %s\n",
                    join( ', ', @differs ) );
            printf( "==> In table:\n%s\n", Dumper( \%row ) );
            printf( "==> In master:\n%s\n",
                    Dumper( $master->{$logic_name_lc} ) );
            $Data::Dumper::Terse  = 1;
            $Data::Dumper::Indent = 0;
            $Data::Dumper::Useqq  = 0;

            push(
              @{ $sql{$dbname} },
              sprintf(
                "-- Updating %s for logic_name '%s'\n"
                  . "UPDATE %s ad,\n"
                  . "       %s a\n"
                  . "SET %s\n"
                  . "WHERE a.logic_name = %s\n"
                  . "AND ad.analysis_id = a.analysis_id;\n",
                join( ', ', @differs ),
                $logic_name_lc,
                $dbh->quote_identifier(
                                  undef, $dbname, 'analysis_description'
                ),
                $dbh->quote_identifier( undef, $dbname, 'analysis' ),
                join(
                  ",\n",
                  map {
                    sprintf( "ad.%s = %s",
                             $_,
                             $dbh->quote( $master->{$logic_name_lc}{$_},
                                          $colinfo->{$_}{'DATA_TYPE'} )
                      )
                    } @differs ),
                $dbh->quote( $logic_name_lc, SQL_VARCHAR ) ),
              sprintf(
                "-- previous value(s) were:\n%s\n",
                join(
                  "\n",
                  map {
                    sprintf( "-- %s = '%s'", $_, $row{$_} )
                    } @differs ) ),
              "\n" );

          } ## end if ( scalar(@differs) ...)
        } ## end else [ if ( !exists( $master->...))]

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
