=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <dev@ensembl.org>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

=head1 NAME

Bio::EnsEMBL::UnconventionalTranscriptAssociation - A class representing
an some sort of unconventional association between a gene and a
transcript.

=head1 SYNOPSIS

  $ex = new Bio::EnsEMBL::UnconventionalTranscriptAssociation( $gene,
    $transcript, $type );

=head1 METHODS

=cut

package Bio::EnsEMBL::UnconventionalTranscriptAssociation;
use vars qw(@ISA);
use strict;

use Bio::EnsEMBL::Utils::Exception qw( warning throw deprecate );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );

=head2 new

  Args [1]   : Bio::EnsEMBL::Gene - the gene which is associated.
  Args [2]   : Bio::EnsEMBL::Transcript - the transcript which is associated.
  Args [3]   : String type - the type of assocation, e.g. "antisense",
               "sense_intronic","sense_overlaping_exonic","chimeric_sense_exonic".
  Example    : $uta = new Bio::EnsEMBL::UnconventionalTranscriptAssociation($gene, $transcript, "antisense")
  Description: create an UnconventionalTranscriptAssociation object.
  Returntype : Bio::EnsEMBL::UnconventionalTranscriptAssociation.
  Exceptions : Wrong argument types
  Caller     : general
  Status     : At risk

=cut

sub new {

  my ($class, $transcript, $gene, $type) = @_;

  $class = ref $class || $class;

  my $self = {};

  if( !ref $gene || ! $gene->isa("Bio::EnsEMBL::Gene") ) {
    throw("$gene is not a Bio::EnsEMBL::Gene!");
  }

  if( !ref $transcript || ! $transcript->isa("Bio::EnsEMBL::Transcript") ) {
    throw("$transcript is not a Bio::EnsEMBL::Transcript!");
  }

  $self->{'gene'} = $gene;
  $self->{'transcript'} = $transcript;
  $self->{'type'} = $type;

  return bless $self, $class;
}


=head2 gene

  Args       : none
  Example    : $gene = $uta->gene()
  Description: Getter/setter for the gene part of this association.
  Returntype : Bio::EnsEMBL::Gene
  Exceptions : none
  Caller     : general
  Status     : At risk

=cut

sub gene {
  my ($self) = shift;

  $self->{'gene'} = shift if (@_);
  return $self->{'gene'};

}

=head2 transcript

  Args       : none
  Example    : $transcript = $uta->transcript()
  Description: Getter/setter for the transcript part of this association.
  Returntype : Bio::EnsEMBL::Transcript
  Exceptions : none
  Caller     : General
  Status     : At risk

=cut

sub transcript {
  my ($self) = shift;

  $self->{'transcript'} = shift if (@_);
  return $self->{'transcript'};

}

=head2 interaction_type

  Args       : none
  Example    : $type = $uta->interaction_type()
  Description: Getter/setter for the interaction_type of this association.
  Returntype : String
  Exceptions : none
  Caller     : General
  Status     : At risk

=cut

sub interaction_type {
  my ($self) = shift;

  $self->{'interaction_type'} = shift if (@_);
  return $self->{'interaction_type'};

}

1;


