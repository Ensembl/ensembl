=head1 LICENSE

See the NOTICE file distributed with this work for additional information
regarding copyright ownership.
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

=head1 NAME

XrefParser::MGI_Desc_Parser

=head1 DESCRIPTION

A parser class to parse the MGI (descriptions) source. Creates 'MISC' xref using MGI accession with description and
also creates the synonyms extracted from the pipe seperated synonym_field

-species = mus_musculus
-species_id = 10090
-data_uri = http://www.informatics.jax.org/downloads/reports/MRK_List2.rpt
-file_format = TSV
-columns = [accession chromosome position start end strand label status marker marker_type feature_type synonym_field]

=head1 SYNOPSIS

  my $parser = XrefParser::MGI_Desc_Parser->new($db->dbh);
  $parser->run({
    source_id  => 58,
    species_id => 10090,
    files      => ["MRK_List2.rpt"],
  });

=cut

package XrefParser::MGI_Desc_Parser;

use strict;
use warnings;
use Carp;
use Text::CSV;

use parent qw( XrefParser::BaseParser );

my $EXPECTED_NUMBER_OF_COLUMNS = 12;



=head2

The run method does the actual parsing and creation of xrefs and synonyms.
Parser gets initialized as noted above and run is called from
Bio::EnsEMBL::Production::Pipeline::Xrefs::ParseSource

=cut

sub run {

  my ($self, $ref_arg) = @_;
  my $source_id  = $ref_arg->{source_id};
  my $species_id = $ref_arg->{species_id};
  my $files      = $ref_arg->{files};
  my $verbose    = $ref_arg->{verbose} // 0;
  my $dbi        = $ref_arg->{dbi} // $self->dbi;

  if ( ( !defined $source_id ) or
       ( !defined $species_id ) or
       ( !defined $files ) ) {
    confess 'Need to pass source_id, species_id and files as pairs';
  }

  my $file = @{$files}[0];

  my $mgi_io = $self->get_filehandle($file);
  if ( !defined $mgi_io ) {
    confess "Could not open $file\n";
  }

  my $input_file = Text::CSV->new({
    sep_char           => "\t",
    quote_char         => undef,
    escape_char        => undef,
    strict             => 1,
  }) or confess "Cannot use file $file: " . Text::CSV->error_diag();

  my $xref_count = 0;
  my $syn_count  = 0;
  my %acc_to_xref;

  # read and validate header
  if ( ! is_file_header_valid( $input_file->header( $mgi_io ) ) ) {
    confess "Malformed or unexpected header in MGI_Desc file '${file}'";
  }

  while ( my $data = $input_file->getline($mgi_io) ) {
    my $accession = $data->[0];
    my $marker = $data->[8];

    $acc_to_xref{$accession} = $self->add_xref({
      acc        => $accession,
      label      => $data->[6],
      desc       => $marker,
      source_id  => $source_id,
      species_id => $species_id,
      dbi        => $dbi,
      info_type  => 'MISC',
    });
    if ( $verbose && !$marker ) {
      print "$accession has no description\n";
    }
    $xref_count += 1;

    if ( defined $acc_to_xref{$accession} ) {
      my @synonyms;
      my $synonym_field = $data->[11];
      if ( $synonym_field ) {
        @synonyms = split qr{ [|] }msx, $synonym_field;
      }
      foreach my $syn (@synonyms) {
        $self->add_synonym( $acc_to_xref{$accession}, $syn, $dbi );
        $syn_count += 1;
      }
    }

  } ## end while ( my $data = $input_file...)

  $mgi_io->eof
    || confess "Error parsing file $file: " . $input_file->error_diag();
  $mgi_io->close();

  if ($verbose) {
    print "$xref_count MGI Description Xrefs added\n";
    print "$syn_count synonyms added\n";
  }

  return 0;    #successful
} ## end sub run


=head2 is_file_header_valid

  Arg [1..N] : list of column names provided by Text::CSV::header()
  Example    : if ( ! is_file_header_valid( $csv->header( $fh ) ) {
                 confess 'Bad header';
               }
  Description: Verifies if the header of a MGI_Desc file follows
               expected syntax.
  Return type: boolean
  Exceptions : none
  Caller     : internal
  Status     : Stable

=cut

sub is_file_header_valid {
  my ( @header ) = @_;

  # Don't bother with parsing column names if their number does not
  # match to begin with
  if ( scalar @header != $EXPECTED_NUMBER_OF_COLUMNS ) {
    return 0;
  }

  my @field_patterns
    = (
        'mgi accession id',
        'chr',
        'cm position',
        'genome coordinate start',
        'genome coordinate end',
        'strand',
        'marker symbol',
        'status',
        'marker name',
        'marker type',
        'feature type',
        'marker synonyms (pipe-separated)',
      );

  my $header_field;
  foreach my $pattern (@field_patterns) {
    $header_field = shift @header;
    if ( $header_field ne $pattern ) {
      return 0;
    }
  }

  # If we have made it this far, all should be in order
  return 1;
}


1;
