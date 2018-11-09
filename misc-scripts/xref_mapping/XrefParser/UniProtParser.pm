=head1 LICENSE

See the NOTICE file distributed with this work for additional
information regarding copyright ownership.

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



package XrefParser::UniProtParser;

use strict;
use warnings;

use Carp;
use Readonly;

use XrefParser::UniProtParser::Extractor;
use XrefParser::UniProtParser::Transformer;
use XrefParser::UniProtParser::Loader;

use parent qw( XrefParser::BaseParser );


Readonly my $DEFAULT_LOADER_BATCH_SIZE => 1000;

Readonly my %source_name_for_section => (
                                         'Swiss-Prot' => 'Uniprot/SWISSPROT',
                                         'TrEMBL'     => 'Uniprot/SPTREMBL',
                                       );


sub run {
  my ( $self, $ref_arg ) = @_;

  my $general_source_id = $ref_arg->{source_id};
  my $species_id        = $ref_arg->{species_id};
  my $species_name      = $ref_arg->{species};
  my $files             = $ref_arg->{files};
  my $release_file      = $ref_arg->{rel_file};
  my $verbose           = $ref_arg->{verbose} // 0;
  my $dbh               = $ref_arg->{dbi} // $self->dbi;

  my $loader_batch_size
    = $ref_arg->{loader_batch_size} // $DEFAULT_LOADER_BATCH_SIZE;

  if ( ( !defined $general_source_id ) or
       ( !defined $species_id ) or
       ( !defined $files ) )
  {
    croak "Need to pass source_id, species_id and files as pairs";
  }

  my $extractor = XrefParser::UniProtParser::Extractor->new({
    'file_names' => $files,
  });
  my $transformer = XrefParser::UniProtParser::Transformer->new({
    'dbh'        => $dbh,
    'species_id' => $species_id,
  });
  my $loader = XrefParser::UniProtParser::Loader->new({
    'batch_size' => $loader_batch_size,
    'dbh'        => $dbh,
  });

  my $source_id_map = $transformer->get_source_id_map();
  if ( $verbose ) {
    while ( my ( $section, $pri_ref ) = each %{ $source_id_map } ) {
      while ( my ( $priority, $source_id ) = each %{ $pri_ref } ) {
        # FIXME: do we really care about the file name here?
        print "$section $priority source id for $files->[0]: $source_id\n";
      }
    }
  }

 RECORD:
  while ( $extractor->get_uniprot_record() ) {

    my $extracted_record = $extractor->extract();

    my $transformed_data
      = $transformer->transform( $extracted_record );

    $loader->load( $transformed_data );

  }

  $loader->flush();
  $extractor->close_input();

  # Extract release numbers from the release file, if provided
  if ( defined $release_file ) {
    my $release_numbers = $self->_get_release_numbers_from_file( $release_file,
                                                                 $verbose );
    $self->_set_release_numbers_on_uniprot_sources( $source_id_map,
                                                    $release_numbers,
                                                    $dbh );
  }

  return 0;
}


# Extract Swiss-Prot and TrEMBL release info from the release file
# Note: the only reason for this to be a method rather than a
# standalone function is that it uses BaseParser::get_filehandle().
sub _get_release_numbers_from_file {
  my ( $self, $release_file_name, $verbose ) = @_;

  my $release_io = $self->get_filehandle( $release_file_name );
  my $release_numbers = {};

  while ( my $line = $release_io->getline() ) {
    my ( $section, $release )
      = ( $line =~ m{
                      \A
                      UniProtKB/
                      (
                        Swiss-Prot
                      |
                        TrEMBL
                      )
                      \s+
                      Release
                      \s+
                      ( [^\n]+ )
                  }msx );
    if ( defined $section ) {
      $release_numbers->{ $source_name_for_section{$section} } = $release;
      if ( $verbose ) {
        print "$section release is '$release'\n";
      }
    }
  }

  $release_io->close();

  return $release_numbers;
}

sub _set_release_numbers_on_uniprot_sources {
  my ( $self, $source_id_map, $release_numbers, $dbh ) = @_;

  foreach my $source ( keys %{ $source_id_map } ) {
    # Priority names are not important here, we only need source IDs
    foreach my $source_id ( values %{ $source_id_map->{$source} } ) {
      $self->set_release( $source_id, $release_numbers->{$source}, $dbh );
    }
  }

  return;
}

1;
