=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut

package XrefParser::DBASSParser;

use strict;
use warnings;

use DBI;
use Carp;

use parent qw( XrefParser::BaseParser );

# This parser will read direct xrefs from a simple comma-delimited file downloaded from the DBASS database.
# The columns of the file should be the following:
#
# 1)    DBASS Gene ID
# 2)    DBASS Gene Name
# 3)    Ensembl Gene ID


sub run {
  my ( $self, $ref_arg ) = @_;
  my $source_id  = $ref_arg->{source_id};
  my $species_id = $ref_arg->{species_id};
  my $files      = $ref_arg->{files};
  my $verbose    = $ref_arg->{verbose};
  my $dbi        = $ref_arg->{dbi};
  $dbi = $self->dbi unless defined $dbi;

  if ( ( !defined $source_id ) or
       ( !defined $species_id ) or
       ( !defined $files ) )
  {
    croak "Need to pass source_id, species_id and files as pairs";
  }
  $verbose |= 0;

  my $filename = @{$files}[0];

  my $file_io = $self->get_filehandle($filename);
  if ( !defined($file_io) ) {
    return 1;
  }

  my $parsed_count = 0;

  # Skip the header line
  $file_io->getline();

  while ( defined( my $line = $file_io->getline() ) ) {

    # strip trailing whitespace
    $line =~ s/\s*$//x;
    # csv format can come with quoted columns, remove them. Note that this
    # assumes actual column values will never contain commas, which while
    # true at present brings into question why bother quoting them in the
    # first place.
    $line =~ s/"//gx;

    my ( $dbass_gene_id, $dbass_gene_name, $ensembl_id, @split_overflow ) =
      split( /,/x, $line );

    if ( !defined($dbass_gene_id) || !defined($ensembl_id) ) {
      printf STDERR ( "Line %d has fewer than three columns.\n",
                      1 + $parsed_count );
      print STDERR ("The parsing failed\n");
      return 1;
    }
    elsif ( scalar @split_overflow > 0 ) {
      printf STDERR ( "Line %d has more than three columns.\n",
                      1 + $parsed_count );
      print STDERR ("The parsing failed\n");
      return 1;
    }

    # Three options here:
    #  1. gene has only one name;
    #  2. synonyms are slash-separated;
    #  3. the second name follows the first one in brackets.
    # FIXME: .* is seriously inefficient because here it results in massive
    # amounts of backtracking. Could we be more specific, i.e. assume
    # some specific format of DBASS names?
    my ( $first_gene_name, $second_gene_name );
    if ( ( $dbass_gene_name =~ m{
                                  (.*)
                                  \s?\/\s?  # typically no ws here but just in case
                                  (.*)
                                }x ) ||
         ( $dbass_gene_name =~ m{
                                  (.*)
                                  \s?  # typically there IS a space before the ( here
                                  [(] (.*) [)]
                                }x ) ) {
      $first_gene_name  = $1;
      $second_gene_name = $2;
    }
    else {
      $first_gene_name = $dbass_gene_name;
      $second_gene_name = undef;
    }

    ++$parsed_count;

    # Do not attempt to create unmapped xrefs
    if ( $ensembl_id ne '' ) {
      my $label       = $first_gene_name;
      my $synonym     = $second_gene_name;
      my $type        = 'gene';
      my $version     = '1';

      my $xref_id =
        $self->get_xref( $dbass_gene_id, $source_id, $species_id, $dbi );

      if ( !defined($xref_id) || $xref_id eq '' ) {
        $xref_id = $self->add_xref(
                                   { acc        => $dbass_gene_id,
                                     version    => $version,
                                     label      => $label,
                                     desc       => undef,
                                     source_id  => $source_id,
                                     dbi        => $dbi,
                                     species_id => $species_id,
                                     info_type  => "DIRECT" } );
      }

      if ( defined($synonym) ) {
        $self->add_synonym( $xref_id, $synonym, $dbi );
      }

      $self->add_direct_xref( $xref_id, $ensembl_id, $type, undef, $dbi );
    }

  } ## end while ( defined( my $line...))

  printf( "%d direct xrefs succesfully parsed\n", $parsed_count )
    if ($verbose);

  $file_io->close();

  return 0;
} ## end sub run


1;
