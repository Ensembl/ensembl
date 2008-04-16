#!/software/bin/perl -w

# $Id$

########################################################################
#                                                                      #
# This script will take the 'xref_config.ini' configuration            #
# file (or whatever file name given on the command line) and           #
# convert it into a SQL file that can be used in place of the old      #
# 'populate_metadata.sql' file found in the 'sql' subdirectory.        #
#                                                                      #
# The output from this script should be redirected to a file that      #
# you manually run to populate your Xref database, just as was done    #
# with 'populate_metadata.sql'.  The safest thing to do is just to     #
# overwrite 'sql/populate_metadata.sql' with the output of this        #
# script.  This will ensure that 'xref_parser.pl populates the Xref    #
# database with the correct data.                                      #
#                                                                      #
########################################################################

use strict;
use warnings;

use Config::IniFiles;

my $config =
  Config::IniFiles->new( -file =>
      ( defined $ARGV[0] && -f $ARGV[0] ? $ARGV[0] : 'xref_config.ini' )
  );

# Do the species.

print( '#' x 80, "\n" );
print("# SPECIES\n");
print("\n");

foreach my $section ( $config->GroupMembers('species') ) {
  my $species_name = substr( $section, 8 );

  my @taxonomy_ids =
    split( /\n/, $config->val( $section, 'taxonomy_id' ) );

  my $species_id = $taxonomy_ids[0];

  printf( "# Species '%s' (id = %d)\n", $species_name, $species_id );

  foreach my $taxonomy_id (@taxonomy_ids) {
    print(   "INSERT INTO species "
           . "(species_id, taxonomy_id, name, aliases)\n" );

    printf( "VALUES (%d, %d, '%s', '%s');\n",
            $species_id, $taxonomy_id, $species_name,
            $config->val( $section, 'aliases' ) );
  }

  print("\n");
}

# Do the sources.

print( '#' x 80, "\n" );
print("# SOURCES\n");
print("\n");

my $source_id = 0;
foreach my $source_section ( sort( $config->GroupMembers('source') ) ) {
  my ( $spaces, $source_name ) =
    $source_section =~ /^source(\s+)(\S+)\s*$/;

  if ( length($spaces) > 1 ) {
    die(
         sprintf( "Too many spaces between the words 'source' and '%s' "
                    . "when reading source section '[%s]'\n",
                  $source_name, $source_section
         ) );
  }

  $config->newval( $source_section, 'id', ++$source_id );

  printf( "# Source '%s' (id = %d)\n", $source_name, $source_id );

  print(   "INSERT INTO source "
         . "(name, source_release, download, ordered, "
         . "priority, priority_description)\n" );

  printf( "VALUES ('%s', '1', '%s', %d, %d, '%s');\n",
          $config->val( $source_section, 'name' ),
          $config->val( $source_section, 'download' ),
          $config->val( $source_section, 'order' ),
          $config->val( $source_section, 'priority' ),
          $config->val( $source_section, 'prio_descr' ) );

  print("\n");
} ## end foreach my $source_section ...

# Do the data files.

print( '#' x 80, "\n" );
print("# DATA FILES\n");
print("\n");

foreach my $species_section ( sort( $config->GroupMembers('species') ) )
{
  my ( $spaces, $species_name ) =
    $species_section =~ /^species(\s+)(\S+)\s*$/;

  if ( length($spaces) > 1 ) {
    die(
         sprintf("Too many spaces between the words 'species' and '%s' "
                   . "when reading species section '[%s]'\n",
                 $species_name, $species_section
         ) );
  }

  my @taxonomy_ids =
    split( /\n/, $config->val( $species_section, 'taxonomy_id' ) );

  my $species_id = $taxonomy_ids[0];

  print( '#', '-' x 79, "\n" );
  printf( "# Data for species '%s' (id = %d)\n",
          $species_name, $species_id );
  print( '#', '-' x 79, "\n" );
  print("\n");

  foreach my $source_name (
     sort( split( /\n/, $config->val( $species_section, 'source' ) ) ) )
  {
    my $source_section = sprintf( "source %s", $source_name );

    if ( !defined( $config->val( $source_section, 'id' ) ) ) {
      die(
           sprintf( "Can not find source section '[%s]'\n"
                      . "while reading species section '[%s]'\n",
                    $source_section, $species_section
           ) );
    }

    printf( "# Data from source '%s' (id = %d)\n",
            $source_name, $config->val( $source_section, 'id' ) );

    print(   "INSERT INTO source_url "
           . "(source_id, species_id, url, release_url, "
           . "file_modified_date, upload_date, parser)\n" );

    my @uris =
      split( /\n/, $config->val( $source_section, 'data_uri', '' ) );

    my $release_uri = $config->val( $source_section, 'release_uri' );

    if ( $release_uri !~ /\w/ ) {
      $release_uri = '\N';
    } else {
      $release_uri = "'$release_uri'";
    }

    printf( "VALUES (%d, %d, '%s', %s, now(), now(), '%s');\n",
            $config->val( $source_section, 'id' ),
            $species_id,
            join( ' ', @uris ),
            $release_uri,
            $config->val( $source_section, 'parser' ) );

    print("\n");
  } ## end foreach my $source_name ( sort...
} ## end foreach my $species_section...
