=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
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

=head1 NAME Bio::EnsEMBL::Intron - A class representing an Intron

=head1 SYNOPSIS

  $ex = new Bio::EnsEMBL::Intron( exon1, exon2 );

=cut


package Bio::EnsEMBL::Intron;
use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::Feature;
use Bio::Seq; # introns have to have sequences...

use Bio::EnsEMBL::Utils::Exception qw( warning throw deprecate );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );


@ISA = qw(Bio::EnsEMBL::Feature);

=head2 new

  Args       : exon1, exon2. The two exons to build the Intron from.
  Example    : $intron = new Bio::EnsEMBL::Intron($exon1, $exon2)
  Description: create an Intron object from two exons.
  Returntype : Bio::EnsEMBL::Intron
  Exceptions : exons not on the same strand or slice.
  Caller     : general
  Status     : Stable

=cut

sub new {
  my ($class,$e1,$e2) = @_;

  $class = ref $class || $class;

  my $self = $class->SUPER::new();

  if($e1->strand == -1){
    $self->{'end'} = ($e1->start)-1;
    $self->{'start'} = ($e2->end)+1;
  }
  else{
    $self->{'start'}= ($e1->end)+1;
    $self->{'end'} = ($e2->start)-1;
  }

  if($e1->strand != $e2->strand){
  #  throw("Exons on different strand. Not allowed");
  }
  else{
    $self->{'strand'} = $e1->strand;
  }

  if($e1->slice ne $e2->slice){
    if($e1->slice->seq_region_name ne $e2->slice->seq_region_name){
      throw("Exons on different slices. Not allowed");
    }
    else{
      warn("Exons have different slice references to the same seq_region\n");
    }
  }
  else{ 
    $self->{'slice'} = $e1->slice;
  }

  $self->{'prev'} = $e1;
  $self->{'next'} = $e2;

  return $self;
}


=head2 prev_Exon

  Args       : none
  Example    : $exon = $intron->prev_Exon
  Description: Returns the exon before this Intron
  Returntype : Bio::EnsEMBL::Exon
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub prev_Exon {
  my ($self) = shift;

  return $self->{'prev'};
}


=head2 next_Exon

  Args       : none
  Example    : $exon = $intron->next_Exon
  Description: Returns the exon after this Intron
  Returntype : Bio::EnsEMBL::Exon
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub next_Exon {
  my ($self) = shift;

  return $self->{'next'};
}


1;


