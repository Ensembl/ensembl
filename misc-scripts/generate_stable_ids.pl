#!/usr/bin/env perl
# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


# Generate stable IDs for genes/transcripts/translations/exons that have none
# Start from current max stable ID + 1

use strict;
use warnings;

use DBI;
use Getopt::Long;

my $help = 0;
my $port = 3306;
my ( $host, $dbname, $user, $pass, @types, $start, $verbose );

if (
  !GetOptions(
    'dbuser|user=s' => \$user,
    'dbpass|pass=s' => \$pass,
    'dbhost|host=s' => \$host,
    'dbport|port=i' => \$port,
    'dbname=s'      => \$dbname,
    'types=s'       => \@types,
    'start=s'       => \$start,
    # USE ENS000001 or ENS for human, ENSMUS00001 or ENSMUS for mouse etc.
    # Don't add G/T/E/P for specific types !!!
    'help!'    => \$help,
    'verbose!' => \$verbose, )
  || (    $help
       || !defined($user)
       || !defined($host)
       || !defined($dbname) ) )
{
  usage();
  exit;
}

if ( !@types ) {
  @types = ( 'gene', 'transcript', 'translation', 'exon' );
}
@types = split( /,/, join( ',', @types ) );

my $dbi =
  DBI->connect( "DBI:mysql:host=$host:port=$port;database=$dbname",
                $user, $pass, { 'RaiseError' => 1 } )
  || die "Can't connect to database\n";

foreach my $type (@types) {
  my $sth;

  # Get starting stable ID, either specified or current max.

  my $new_stable_id;
  if ( defined($start) ) {
    $new_stable_id = $start;
  } else {
    $new_stable_id = get_highest_stable_id( $dbi, $type );
  }

  if ($verbose) {
    print("Highest, pruned $type\_stable_id found : $new_stable_id \n");
  }

  # Get timestamp so all new stable IDs have the same created/modified
  # dates.
  $sth = $dbi->prepare("SELECT NOW()");
  $sth->execute();
  my $ts;
  if ( my @row = $sth->fetchrow_array() ) {
    $ts = $row[0];
  } else {
    die("Can't get timestamp\n");
  }
  $sth->finish();

  # get a list of objects that don't currently have stable IDs assigned
  # and assign new ones, incrementing & preserving formatting as we go
  my $sql =
      "SELECT ${type}_id "
    . "FROM $type "
    . "WHERE stable_id IS NULL";
  $sth = $dbi->prepare($sql);
  $sth->execute();

  while ( my @row = $sth->fetchrow_array() ) {
    ( $new_stable_id, my $nis ) =
      @{ increment_stable_id( $new_stable_id, $type ) };
    print(   "UPDATE $type SET stable_id = \'$nis\', version = 1, created_date = \'$ts\', modified_date = \'$ts\'"
           . " WHERE ${type}_id = $row[0];\n" );
  }
} ## end foreach my $type (@types)

#-------------------------------------------------------------------------------

sub increment_stable_id {
  my ( $stable_id, $type ) = @_;

  my ( $prefix, $suffix );

  # Check stable_id format ...
  if ( $stable_id =~ m/([a-zA-Z]+)([0-9]+)/ ) {
    ( $prefix, $suffix ) = $stable_id =~ /([a-zA-Z]+)([0-9]+)/;
  } elsif ( $stable_id =~ m/([a-zA-Z]+)/ ) {
    $prefix = $stable_id;
  } else {
    die(   "unrecongnized stable_id format "
         . "- should match ([a-zA-Z]+)([0-9]+) or ([a-zA-Z]+) !!\n" );
  }

  my $new_sid;
  if ( $type eq 'gene' ) {
    $new_sid = $prefix . 'G';
  } elsif ( $type eq 'transcript' ) {
    $new_sid = $prefix . 'T';
  } elsif ( $type eq 'translation' ) {
    $new_sid = $prefix . 'P';
  } elsif ( $type eq 'exon' ) {
    $new_sid = $prefix . 'E';
  }
  my $new_stable_id = sprintf( "%s%011d", $new_sid, $suffix + 1 );

  my $old = sprintf( "%s%011d", $prefix, $suffix + 1 );

  return [ $old, $new_stable_id ];
} ## end sub increment_stable_id

#-------------------------------------------------------------------------------

sub get_max_stable_id_from_gene_archive {
  my ( $dbi, $type ) = @_;

  # Try to get from relevant archive.
  my $sth = $dbi->prepare("SELECT MAX($type) FROM gene_archive");
  $sth->execute();

  my $rs;
  if ( my @row = $sth->fetchrow_array ) {
    $rs = $row[0];
  }

  if ( length($rs) <= 0 ) {
    print( STDERR "no entry for $type found in gene_archive table "
           . "- returning undef\n" );
    return undef;
  }

  return $rs;
}

#-------------------------------------------------------------------------------

sub get_highest_stable_id {
  my ( $dbi, $type ) = @_;

  my ( $highest_from_current, $highest_from_archive );

  # Get highest stable ID from the relevant table.

  my $sth = $dbi->prepare("SELECT MAX(stable_id) FROM $type");
  $sth->execute();

  if ( my @row = $sth->fetchrow_array() ) {
    $highest_from_current = $row[0];
  } else {
    die("Can't get max $type stable ID from $type\n");
  }

  if ( length($highest_from_current) == 0 ) {
    print( STDERR
           " Warning ! length of stable_id for $type is zero \n" );
  }

  if ( $type eq "exon" ) {
    # Archive doesn't store information about exon_stable_ids so try
    # without archive first.

    if ( length($highest_from_current) == 0 ) {
      print( "\n"
        . "WARNING:\n"
        . "No stable_id for exon found \n"
        . "I got no prefix to generate new stable_ids for type $type!!! "
        . "- I'll try to use gene_archive now\n" );

      my $max =
        get_max_stable_id_from_gene_archive( $dbi, "gene_stable_id" );

      my $prefix;

      if ( length($max) > 0 ) {
        ( $prefix, my $suffix ) = $max =~ /([a-zA-Z]+)([0-9]+)/;
        $prefix =~ s/G$//g;
      } else {
        die(  "ERROR: "
            . "No entries in table exon and "
            . "gene_archive tables found\n"
            . "Don't know which species prefix to use for species.\n" );

        $highest_from_current = sprintf( "%s%011d", $prefix, 0 );
      }
    } ## end if ( length($highest_from_current...))

    return $highest_from_current;
  } ## end if ( $type eq "exon" )

  # and from relevant archive

  $highest_from_archive =
    get_max_stable_id_from_gene_archive( $dbi, $type . "_stable_id" );

  my $max =
    ( $highest_from_current ge $highest_from_archive )
    ? $highest_from_current
    : $highest_from_archive;

  if ( length($max) == 0 ) {
    die(   "ERROR: "
         . "No stable_id in table gene_archive "
         . "or found in $type\_stable_id - tables\n" );
  }

  # Assuming that this is a correctly formatted stable id -> remove the
  # G/T/P/E for exon etc.

  my ( $prefix, $suffix ) = $max =~ /([a-zA-Z]+)([0-9]+)/;
  if ( $type eq 'exon' ) {
    $prefix =~ s/E$//;
  } elsif ( $type eq 'gene' ) {
    $prefix =~ s/G$//;
  } elsif ( $type eq 'transcript' ) {
    $prefix =~ s/T$//;
  } elsif ( $type eq 'translation' ) {
    $prefix =~ s/P$//;
  }

  return $prefix . $suffix;
} ## end sub get_highest_stable_id

#-------------------------------------------------------------------------------

sub usage {
  print <<USAGE_END;

  USAGE:

  generate_stable_ids.pl -dbuser|user {user}
                         -dbpass|pass {password}
                         -dbhost|host {host}
                         -dbport|port {port}
                         -dbname {database}
                         -types {gene,exon,transcript,translation}
                         -start {first stable ID}

  Argument to -types is a comma-separated list of types of stable IDs to
  be produced.

  If the -types argument is ommitted, stable IDs are generated for all
  types (gene,transcript,translation,exon).

  Assigns stable IDs to objects that currently have none.  The starting
  stable ID is found by incrementing the highest current stable ID for
  that type *or* by using -start argument.  If no -start option is used
  the script tries to find the latest given stable_id for each object
  by looking up the <OBJ>_stable_id tables in the database and the
  gene_archive table (only for gene, translation and transcript, not for
  exon stable IDs!)

  Note:

  The -start option requires to not submit an initial stable ID without
  any gene, transcript, translation or exon specifier, as in

   ENSMUS000001 (not ENSMUSG0001, then you end up with stable IDs like
                ENSMUSGG001, ENSMUSGT0001.)

  Again, the parameter to -start should be the stable ID you wish
  to start from without the gene, transcript, translation or exon
  specifier.

  Examples:

    - to generate only exon stable IDs starting with 223 for mouse, use

                    -start ENSMUS222 -types exon

    - to generate exon and gene stable IDs, starting with 223 for mouse,
      use

                    -start ENSMUS222 -types exon,gene

    - to generate all types of stable IDs for human, which all start
      with ID 666 use

                    -start ENS665

    - to generate a whole new set of stable_ids (exon, transcript,
      translation, gene) starting with 1 for an organism with prefix
      ENSNEW you can use one of the following options :

                    -start ENSNEW0          <or>
                    -start ENSNEW0          <or>
                    -start ENSNEW00000000

  The script produces SQL which can be run against the target database.

USAGE_END
} ## end sub usage
