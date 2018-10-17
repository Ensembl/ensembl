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

package XrefParser::HPAParser;

use strict;
use warnings;
use Carp;
use base qw( XrefParser::BaseParser);

# This parser will read direct xrefs from a simple comma-delimited file downloaded from the Human Protein Atlas (HPA) database.
# The database contains two types of antibody, their own HPA antibodies and Collaborator antibody (CAB) commercial antibodies. 
# The columns of the file should be the following:
#
# 1)    Antibody
# 2)    Antibody ID
# 3)    Ensembl Peptide ID
# 4)	Link (URL)


sub run {
  my ($self, $ref_arg) = @_;
  my $source_id    = $ref_arg->{source_id};
  my $species_id   = $ref_arg->{species_id};
  my $files        = $ref_arg->{files};
  my $verbose      = $ref_arg->{verbose};
  my $dbi          = $ref_arg->{dbi};
  $dbi = $self->dbi unless defined $dbi;

  croak "Need to pass source_id, species_id, files and rel_file as pairs"
    unless defined $source_id and defined $species_id and defined $files;

  $verbose |= 0;
  my $file = @{$files}[0];

  my $file_io = $self->get_filehandle($file);
  croak "ERROR: Could not open $file\n" unless defined $file_io;

  my $parsed_count = 0;
  $file_io->getline(); # skip header # Antibody,antibody_id,ensembl_peptide_id,link

  while ( defined( my $line = $file_io->getline() ) ) {
    $line =~ s/\s*$//;
		
    my ( $antibody, $antibody_id, $ensembl_peptide_id, $link ) = split( /,/, $line );

    croak sprintf("Line %d contains has less than two columns.\nParsing failed", 1 + $parsed_count)
      unless defined $antibody and defined $ensembl_peptide_id;
	
    my $label       = $antibody;
    my $type        = 'translation';
    my $version     = '1';

    ++$parsed_count;

    my $xref_id = $self->get_xref( $antibody_id, $source_id, $species_id, $dbi );
    $xref_id = $self->add_xref({ acc        => $antibody_id,
				 version    => $version,
				 label      => $label,
				 source_id  => $source_id,
				 species_id => $species_id,
				 dbi        => $dbi,
				 info_type  => "DIRECT"} ) unless defined $xref_id or $xref_id ne '';

    $self->add_direct_xref( $xref_id, $ensembl_peptide_id, $type, undef, $dbi);
  } ## end while ( defined( my $line...

  printf( "%d direct xrefs succesfully parsed\n", $parsed_count ) if $verbose;

  $file_io->close();

  return 0;
} ## end sub run

1;
