
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

package XrefParser::MGI_Desc_Parser;

use strict;
use warnings;
use Carp;
use File::Basename;
use Text::CSV;

use parent qw( XrefParser::BaseParser );

sub run {

  my ( $self, $ref_arg ) = @_;
  my $source_id  = $ref_arg->{source_id};
  my $species_id = $ref_arg->{species_id};
  my $files      = $ref_arg->{files};
  my $verbose    = $ref_arg->{verbose} // 0;
  my $dbi        = $ref_arg->{dbi} // $self->dbi;

  if ( ( !defined $source_id )
    or ( !defined $species_id )
    or ( !defined $files ) )
  {
    croak "Need to pass source_id, species_id and files as pairs";
  }

  my $file = shift @{$files};

  my $mgi_io = $self->get_filehandle($file);

  if ( !defined $mgi_io ) {
    croak "ERROR: Could not open $file\n";
  }

  my $xref_count = 0;
  my $syn_count  = 0;
  my %acc_to_xref;

  my $input_file = Text::CSV->new(
    {
      sep_char           => "\t",
      empty_is_undef     => 1,
      strict             => 1,
      allow_loose_quotes => 1,
    }
  ) or croak "Cannot use file $file: " . Text::CSV->error_diag();

  # expected columns
  my @expected_columns =
    qw(accession chromosome position start end strand label status marker marker_type feature_type synonym_field);

  # read header
  my $header = $input_file->getline($mgi_io);
  if ( scalar @{$header} != scalar @expected_columns ) {
    croak "input file $file has an incorrect number of columns";
  }
  $input_file->column_names( \@expected_columns );

  while ( my $data = $input_file->getline_hr($mgi_io) ) {
    my $accession = $data->{'accession'};
    my $marker = defined( $data->{'marker'} ) ? $data->{'marker'} : undef;
    $acc_to_xref{$accession} = $self->add_xref(
      {
        acc        => $accession,
        label      => $data->{'label'},
        desc       => $marker,
        source_id  => $source_id,
        species_id => $species_id,
        dbi        => $dbi,
        info_type  => "MISC"
      }
    );
    if ( $verbose && !$marker ) {
      print "$accession has no description\n";
    }
    $xref_count++;
    my @synonyms;
    if ( defined( $acc_to_xref{$accession} ) ) {
      @synonyms = split( /\|/, $data->{'synonym_field'} )
        if ( $data->{'synonym_field'} );
      foreach my $syn (@synonyms) {
        $self->add_synonym( $acc_to_xref{$accession}, $syn, $dbi );
        $syn_count++;
      }
    }
  }
  $mgi_io->eof
    or croak "Error parsing file $file: " . $input_file->error_diag();
  $mgi_io->close();

  if ($verbose) {
    print "$xref_count MGI Description Xrefs added\n";
    print "$syn_count synonyms added\n";
  }

  return 0;    #successful
}

1;

