=head1 LICENSE

  Copyright (c) 1999-2010 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 METHODS

=cut

package Bio::EnsEMBL::IdMapping::InternalIdMapper::EnsemblExonGeneric;

use strict;
use warnings;
no warnings 'uninitialized';

use Bio::EnsEMBL::IdMapping::InternalIdMapper::BaseMapper;
our @ISA = qw(Bio::EnsEMBL::IdMapping::InternalIdMapper::BaseMapper);

use Bio::EnsEMBL::Utils::Exception qw(throw warning);

  
#
# basic mapping
#
sub init_basic {
  my $self = shift;
  my $num = shift;
  my $esb = shift;
  my $mappings = shift;
  my $exon_scores = shift;

  $self->logger->info("Basic exon mapping...\n", 0, 'stamped');

  $mappings = $self->basic_mapping($exon_scores, "exon_mappings$num");
  $num++;
  my $new_scores = $esb->create_shrinked_matrix($exon_scores, $mappings,
    "exon_matrix$num");

  return ($new_scores, $mappings);
}


#
# reduce score for mappings of exons which do not belong to mapped
# transcripts
#
sub mapped_transcript {
  my $self = shift;
  my $num = shift;
  my $esb = shift;
  my $mappings = shift;
  my $exon_scores = shift;

  $self->logger->info("Exons in mapped transcript...\n", 0, 'stamped');

  unless ($exon_scores->loaded) {
    $esb->non_mapped_transcript_rescore($exon_scores, $mappings);
    $exon_scores->write_to_file;
  }

  $mappings = $self->basic_mapping($exon_scores, "exon_mappings$num");
  $num++;
  my $new_scores = $esb->create_shrinked_matrix($exon_scores, $mappings,
    "exon_matrix$num");

  return ($new_scores, $mappings);
}
  

1;

