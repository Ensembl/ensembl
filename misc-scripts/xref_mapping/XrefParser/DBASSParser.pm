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

# For non-destructive substitutions in regexps (/r flag)
require 5.014_000;

use Carp;
use List::Util;
use Readonly;
use Text::CSV;

use parent qw( XrefParser::BaseParser );

# This parser will read direct xrefs from a comma-delimited file downloaded from the DBASS Web site.
# The columns of the file should be the following:
#
# 1)    DBASS Gene ID
# 2)    DBASS Gene Name
# 3)    Ensembl Gene ID
#
# where 2) can be either a single name, a 'name/synonym' pair, or a 'name (synonym)' pair.
# Column values, including empty strings, can be surrounded by pairs of double quotes.


Readonly my $EXPECTED_NUMBER_OF_COLUMNS => 3;


sub run {
  my ( $self, $ref_arg ) = @_;
  my $source_id  = $ref_arg->{source_id};
  my $species_id = $ref_arg->{species_id};
  my $files      = $ref_arg->{files};
  my $verbose    = $ref_arg->{verbose} // 0;
  my $dbi        = $ref_arg->{dbi} // $self->dbi;

  if ( ( !defined $source_id ) or
       ( !defined $species_id ) or
       ( !defined $files ) )
  {
    croak 'Need to pass source_id, species_id and files as pairs';
  }
  my $csv = Text::CSV->new()
    || croak 'Failed to initialise CSV parser: ' . Text::CSV->error_diag();

  my $filename = @{$files}[0];

  my $file_io = $self->get_filehandle($filename);
  if ( !defined($file_io) ) {
    croak "Failed to acquire a file handle for '${filename}'";
  }

  if ( ! is_file_header_valid( $csv->getline( $file_io ) ) ) {
    croak "Malformed or unexpected header in DBASS file '${filename}'";
  }

  my $processed_count = 0;
  my $unmapped_count = 0;

  while ( defined( my $line = $csv->getline ( $file_io ) ) ) {

    if ( scalar @{ $line } != $EXPECTED_NUMBER_OF_COLUMNS ) {
      croak 'Line ' . (2 + $processed_count + $unmapped_count)
        . " of input file '${filename}' has an incorrect number of columns";
    }

    # Do not modify the contents of @{$line}, only the output - hence the /r.
    my ( $dbass_gene_id, $dbass_gene_name, $ensembl_id )
      = map { s{\s+\z}{}rmsx } @{ $line };

    # Do not attempt to create unmapped xrefs. Checking truthiness is good
    # enough here because the only non-empty string evaluating as false is
    # not a valid Ensembl stable ID.
    if ( $ensembl_id ) {

      # DBASS files list synonyms in two ways: either "FOO (BAR)" (with or
      # without space) or "FOO/BAR". Both forms are relevant to us.
      my ( $first_gene_name, $second_gene_name );
      if ( ( $dbass_gene_name =~ m{
                                    (.*)
                                    \s?\/\s?  # typically no ws here but just in case
                                    (.*)
                                  }msx ) ||
           ( $dbass_gene_name =~ m{
                                    (.*)
                                    \s?  # there are entries both with and without ws
                                    [(] (.*) [)]
                                  }msx ) ) {
        $first_gene_name  = $1;
        $second_gene_name = $2;
      }
      else {
        $first_gene_name = $dbass_gene_name;
        $second_gene_name = undef;
      }

      my $label       = $first_gene_name;
      my $synonym     = $second_gene_name;
      my $type        = 'gene';
      my $version     = '1';

      my $xref_id =
        $self->get_xref( $dbass_gene_id, $source_id, $species_id, $dbi );

      if ( !defined($xref_id) || $xref_id eq q{} ) {
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

      ++$processed_count;
    }
    else {
      ++$unmapped_count;
    }

  } ## end while ( defined( my $line...))

  $csv->eof || croak 'Error parsing CSV: ' . $csv->error_diag();
  $file_io->close();

  if ($verbose) {
    printf( "%d direct xrefs succesfully processed\n", $processed_count );
    printf( "Skipped %d unmapped xrefs\n", $unmapped_count );
  }

  return 0;
} ## end sub run


=head2 is_file_header_valid

  Arg [1]    : String file header line
  Example    : if (!is_file_header_valid($header_line)) {
                 croak 'Bad header';
               }
  Description: Verifies if the header of a DBASS file follows expected
               syntax and contains expected column names.
  Return type: boolean
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub is_file_header_valid {
  my ( $header ) = @_;

  # Don't bother with parsing column names if their number does not
  # match to begin with
  if ( scalar @{ $header } != $EXPECTED_NUMBER_OF_COLUMNS ) {
    return 0;
  }

  my @fields_ok;

  my ( $dbass_end ) = ( $header->[0] =~ m{ DBASS (3|5) GeneID }msx );
  push @fields_ok, defined $dbass_end;

  my $dbass_name_ok = ( $header->[1] =~ m{ DBASS ${dbass_end} GeneName }msx );
  push @fields_ok, $dbass_name_ok;

  my $ensembl_id_ok = ( $header->[2] eq 'EnsemblGeneNumber' );
  push @fields_ok, $ensembl_id_ok;

  # All fields must be in order
  return List::Util::all { $_ } @fields_ok;
}


1;
