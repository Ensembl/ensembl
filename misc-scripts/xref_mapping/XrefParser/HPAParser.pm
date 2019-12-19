=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

XrefParser::HPAParser

=head1 DESCRIPTION

This parser will read and creates direct xrefs from a simple comma-delimited file downloaded from the Human Protein Atlas (HPA) database.
The database contains two types of antibody, their own HPA antibodies and Collaborator antibody (CAB) commercial antibodies.

 data_uri        = http://www.proteinatlas.org/download/xref.php

The columns of the file should be the following:

 1)  Antibody
 2)  Antibody ID
 3)  Ensembl Peptide ID
 4)  Link (URL)

 Antibody,antibody_id,ensembl_peptide_id,link
 CAB000001,1,ENSP00000363822,http://www.proteinatlas.org/ENSG00000169083-AR
 CAB000001,1,ENSP00000379358,http://www.proteinatlas.org/ENSG00000169083-AR

=head1 SYNOPSIS

  my $parser = XrefParser::HPAParser->new($db->dbh);
  $parser->run({
    source_id  => 11,
    species_id => 9606,
    files      => ["hpa.txt"],
  });

=cut

package XrefParser::HPAParser;

use strict;
use warnings;

use Carp;
use Text::CSV;

use parent qw( XrefParser::BaseParser);

my $EXPECTED_NUMBER_OF_COLUMNS = 4;



=head2
The run method does the actual parsing and creation of direct xrefs.
Parser gets initialized as noted above and run is called from
Bio::EnsEMBL::Production::Pipeline::Xrefs::ParseSource

my $parser = XrefParser::HPAParser->new($db->dbh);
$parser->run(...);

=cut

sub run {
  my ( $self, $ref_arg ) = @_;
  my $source_id  = $ref_arg->{source_id};
  my $species_id = $ref_arg->{species_id};
  my $files      = $ref_arg->{files};
  my $verbose    = $ref_arg->{verbose} // 0;
  my $dbi        = $ref_arg->{dbi} // $self->dbi;

  if ( ( !defined $source_id ) or
       ( !defined $species_id ) or
       ( !defined $files ) ) {
    confess 'Need to pass source_id, species_id, and files';
  }

  my $file = @{$files}[0];

  my $file_io = $self->get_filehandle($file);
  if ( !defined $file_io ) {
    confess "Could not open $file\n";
  }

  my $input_file = Text::CSV->new({
    sep_char       => q{,},
    empty_is_undef => 1,
    strict         => 1,
  }) or confess "Cannot use file $file: " . Text::CSV->error_diag();

  if ( ! is_file_header_valid( $input_file->header( $file_io ) ) ) {
    confess "Malformed or unexpected header in HPA file '${file}'";
  }

  my $parsed_count = 0;
  while ( my $data = $input_file->getline($file_io) ) {
    my ( $antibody_name, $antibody_id, $ensembl_id ) = @{ $data };

    $self->add_to_direct_xrefs({
      acc        => $antibody_id,
      version    => '1',
      label      => $antibody_name,
      stable_id  => $ensembl_id,
      type       => 'translation',
      source_id  => $source_id,
      species_id => $species_id,
      info_type  => 'DIRECT'
    });

    ++$parsed_count;
  } ## end while

  $input_file->eof or
    confess "Error parsing file $file: " . $input_file->error_diag();
  $file_io->close();

  if ($verbose) {
    printf( "%d direct xrefs succesfully parsed\n", $parsed_count );
  }

  return 0;
} ## end sub run


=head2 is_file_header_valid

  Arg [1..N] : list of column names provided by Text::CSV::getline()
  Example    : if ( ! is_file_header_valid( $csv->getline( $fh ) ) {
                 confess 'Bad header';
               }
  Description: Verifies if the header of a HPA file follows expected
               syntax.
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
        qr{ antibody }msx,
        qr{ antibody_id }msx,
        qr{ ensembl_peptide_id }msx,
        qr{ link }msx,
      );

  my $header_field;
  foreach my $pattern (@field_patterns) {
    $header_field = shift @header;
    # Make sure we run the regex match in scalar context
    return 0 unless scalar ( $header_field =~ m{ $pattern }msx );
  }

  # If we have made it this far, all should be in order
  return 1;
}


1;
