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


# FIXME: this belongs in BaseParser
my $ERR_SOURCE_ID_NOT_FOUND = -1;



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
  if ( $wiki_source_id == $ERR_SOURCE_ID_NOT_FOUND ) {
    confess 'Failed to retrieve WikiGene source ID';
  }

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
  $input_file->column_names( @{ $input_file->getline($eg_io) } );

  my $xref_count = 0;
  my $syn_count  = 0;
  my %seen;    # record already processed xrefs

  # read data and load xrefs
 RECORD:
  while ( my $data = $input_file->getline_hr($eg_io) ) {
    # species_id corresponds to the species taxonomy id, see:
    # https://github.com/Ensembl/ensembl-xref/pull/31#issuecomment-445838474
    if ( $data->{'#tax_id'} ne $species_id ) {
      next RECORD;
    }

    my $acc = $data->{'GeneID'};
    if ( exists $seen{$acc} ) {
      next RECORD;
    }

    my $symbol = $data->{'Symbol'};
    my $desc   = $data->{'description'};

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

    my ( @syn ) = split qr{ \| }msx, $data->{'Synonyms'};
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


1;
