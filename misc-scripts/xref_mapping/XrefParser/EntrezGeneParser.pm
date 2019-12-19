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

XrefParser::EntrezGeneParser

=head1 DESCRIPTION

This parser will read and create dependent xrefs from a simple
comma-delimited file downloaded from the EntrezGene database.

=head1 SYNOPSIS

  my $parser = XrefParser::EntrezGeneParser->new($db->dbh);
  $parser->run({
    source_id  => 11,
    species_id => 9606,
    files      => [ "gene_info.gz" ],
  });

=cut

package XrefParser::EntrezGeneParser;

use strict;
use warnings;

use Carp;
use Text::CSV;

use parent qw( XrefParser::BaseParser );


my $EXPECTED_NUMBER_OF_COLUMNS = 16;



=head2 run

  Arg [1]    : HashRef standard list of arguments from ParseSource
  Description: Add dependent xrefs from EntrezGene to the xref database
  Return type: Int; 0 upon success
  Exceptions : throws on all processing errors
  Caller     : ParseSource in the xref pipeline

=cut

sub run {

  my ( $self, $ref_arg ) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $species_name = $ref_arg->{species};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose} // 0;
  my $dbi          = $ref_arg->{dbi} // $self->dbi;

  if ( ( !defined $source_id ) or
       ( !defined $species_id ) or
       ( !defined $files ) )
    {
    confess 'Need to pass source_id, species_id and files';
  }

  my $file = @{$files}[0];

  my $wiki_source_id =
    $self->get_source_id_for_source_name( 'WikiGene', undef, $dbi );

  my $eg_io = $self->get_filehandle($file);
  if ( !defined $eg_io ) {
    confess "Could not open $file";
  }

  my $input_file = Text::CSV->new({
    sep_char => "\t",
    empty_is_undef => 1,
    allow_loose_quotes => 1
  })
    || confess "Cannot use file $file: " . Text::CSV->error_diag();

  # process header
  if ( ! is_file_header_valid( $input_file->header( $eg_io ) ) ) {
    confess "Malformed or unexpected header in EntrezGene file '${file}'";
  }

  my $xref_count = 0;
  my $syn_count  = 0;
  my %seen;    # record already processed xrefs

  # read data and load xrefs
 RECORD:
  while ( my $data = $input_file->getline($eg_io) ) {
    my ( $tax_id, $acc, $symbol, undef, $synonyms, undef, undef, undef, $desc ) = @{ $data };

    # species_id corresponds to the species taxonomy id, see:
    # https://github.com/Ensembl/ensembl-xref/pull/31#issuecomment-445838474
    if ( $tax_id ne $species_id ) {
      next RECORD;
    }

    if ( exists $seen{$acc} ) {
      next RECORD;
    }

    $self->add_xref({
      acc        => $acc,
      label      => $symbol,
      desc       => $desc,
      source_id  => $source_id,
      species_id => $species_id,
      dbi        => $dbi,
      info_type  => 'DEPENDENT'
    });
    $self->add_xref({
      acc        => $acc,
      label      => $symbol,
      desc       => $desc,
      source_id  => $wiki_source_id,
      species_id => $species_id,
      dbi        => $dbi,
      info_type  => 'DEPENDENT'
    });
    $xref_count += 1;

    my @syn = split qr{ \| }msx, $synonyms;
    foreach my $synonym ( @syn ) {
      if ( $synonym ne q{-} ) {
        $self->add_to_syn( $acc, $source_id, $synonym, $species_id, $dbi );
        $syn_count += 1;
      }
    }

    $seen{$acc} = 1;
  } ## end while ( my $data = $input_file...)

  $input_file->eof ||
    confess "Error parsing file $file, should be EOF: " . $input_file->error_diag();
  $eg_io->close();

  if ( $verbose ) {
    print $xref_count . " EntrezGene Xrefs added with $syn_count synonyms\n";
  }

  return 0;
} ## end sub run


=head2 is_file_header_valid

  Arg [1..N] : list of column names provided by Text::CSV::getline()
  Example    : if ( ! is_file_header_valid( $csv->getline( $fh ) ) {
                 confess 'Bad header';
               }
  Description: Verifies if the header of a EntrezGene file follows expected
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
        qr{ \A [#]? \s* tax_id }msx,
        qr{ geneid }msx,
        qr{ symbol }msx,
        qr{ locustag }msx,
        qr{ synonyms }msx,
        qr{ dbxrefs }msx,
        qr{ chromosome }msx,
        qr{ map_location }msx,
        qr{ description }msx,
        qr{ type_of_gene }msx,
        qr{ symbol_from_nomenclature_authority }msx,
        qr{ full_name_from_nomenclature_authority }msx,
        qr{ nomenclature_status }msx,
        qr{ other_designations }msx,
        qr{ modification_date }msx,
        qr{ feature_type }msx,
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
