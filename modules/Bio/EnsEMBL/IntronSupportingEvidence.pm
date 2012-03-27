package Bio::EnsEMBL::IntronSupportingEvidence;

=pod

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

=head1 NAME

Bio::EnsEMBL::IntronSupportingEvidence

=head1 DESCRIPTION

Used to represent evidence used to delcare the Intron

=head1 METHODS

=cut


use strict;
use warnings;
use base qw/Bio::EnsEMBL::Storable/;

use Bio::EnsEMBL::Utils::Argument qw/rearrange/;
use Bio::EnsEMBL::Utils::Scalar qw/assert_ref/;

=head2 new

  Arg [-ADAPTOR]      : Bio::EnsEMBL::DBSQL::IntronSupportingEvidenceAdaptor
  Arg [-DBID]         : Integer $dbID
  Arg [-INTRON]       : Bio::EnsEMBL::Intron $intron
  Arg [-HIT_NAME]     : String The name of the hit
  Arg [-SCORE]        : Double The score associated with the supporting evidence
  Arg [-SCORE_TYPE]   : String The type of score we are representing
  Example           : Bio::EnsEMBL::IntronSupportingEvidence->new();
  Description       : Returns a new instance of this object
  Returntype        : Bio::EnsEMBL::IntronSupportEvidence
  Exceptions        : Thrown if data is not as expected

=cut

sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args);
  my ($intron, $hit_name, $score, $score_type) = 
    rearrange([qw/intron hit_name score score_type/], @args);
  
  $self->intron($intron);
  $self->hit_name($hit_name);
  $self->score($score);
  $self->score_type($score_type);
  
  return $self;
}

sub intron {
  my ($self, $intron) = @_;
  if(defined $intron) {
    assert_ref($intron, 'Bio::EnsEMBL::Intron', 'intron');
  	$self->{'intron'} = $intron;
  }
  return $self->{'intron'};
}

sub hit_name {
  my ($self, $hit_name) = @_;
  $self->{'hit_name'} = $hit_name if defined $hit_name;
  return $self->{'hit_name'};
}

sub score {
  my ($self, $score) = @_;
  $self->{'score'} = $score if defined $score;
  return $self->{'score'};
}

sub score_type {
  my ($self, $score_type) = @_;
  $self->{'score_type'} = $score_type if defined $score_type;
  return $self->{'score_type'};
}

1;
