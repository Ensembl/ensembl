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

=head1 NAME Bio::EnsEMBL::Intron - A class representing an Intron

=head1 SYNOPSIS

  $intron = Bio::EnsEMBL::Intron->new( exon1, exon2, $analysis );

=cut


package Bio::EnsEMBL::Intron;
use strict;
use warnings;

use Bio::EnsEMBL::Utils::Exception qw( warning throw );

use base qw(Bio::EnsEMBL::Feature);

=head2 new

  Arg [1]    : Bio::EnsEMBL::Exon The 5' exon for the intron; required
  Arg [2]    : Bio::EnsEMBL::Exon The 3' exon for the intron; required
  Arg [3]    : Bio::EnsEMBL::Analysis Analysis to link to this Intron
  Example    : $intron = new Bio::EnsEMBL::Intron($exon1, $exon2)
  Description: Create an Intron object from two exons and an optional analysis
  Returntype : Bio::EnsEMBL::Intron
  Exceptions : exons not on the same strand or slice.
  Caller     : general
  Status     : Stable

=cut

sub new {
  my ( $proto, $e1, $e2, $analysis ) = @_;

  my $class = ref $proto || $proto;

  my $self = $class->SUPER::new();

  if ( $e1->strand() == -1 ) {
    $self->{'end'}   = $e1->start() - 1;
    $self->{'start'} = $e2->end() + 1;
  } else {
    $self->{'start'} = $e1->end() + 1;
    $self->{'end'}   = $e2->start() - 1;
  }

  if ( $e1->strand() != $e2->strand() ) {
    #  throw("Exons on different strand. Not allowed");
  } else {
    $self->{'strand'} = $e1->strand();
  }

  if ( $e1->slice() != $e2->slice() ) {
    if ( ( $e1->slice()->seq_region_name() ne
           $e2->slice()->seq_region_name() )
         && ( $e1->slice()->coord_system_name() ne
              $e2->slice()->coord_system_name() ) )
    {
      throw("Exons on different slices. Not allowed");
    } else {
      warning("Exons have different slice references to the same seq_region");
    }
  } else {
    $self->{'slice'} = $e1->slice();
  }
  
  if($analysis) {
    $self->analysis($analysis);
  }

  $self->{'prev'} = $e1;
  $self->{'next'} = $e2;

  return $self;
} ## end sub new

=head2 length

  Args       : none
  Example    : $length = $intron->length();
  Description: Returns the length of this intron
  Returntype : Integer
  Exceptions : none
  Caller     : general
  Status     : Stable

=cut

sub length {
  my ($self) = @_;

  # TODO: Introns on circular slices, see Feature.pm but allow for
  # zero-length introns.

  return $self->{'end'} - $self->{'start'} + 1;
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

=head2 is_splice_canonical

  Example     : my $canonical = $intron->is_splice_canonical(); 
  Description : Indicates if the splice site is considered normal. This means
                splice site variants equal to (D == donor, A == acceptor)
                  GT (D) => AG (A) 
                  AT (D) => AC (A)
                  GC (D) => AG (A)
  Returntype  : Boolean indicating if the splice was as expected
  Exceptions  : See splice_seq

=cut

sub is_splice_canonical {
  my ($self) = @_;
  my $splice = join q{}, @{$self->splice_seq()}; 
  my $canonical = {
    'GTAG' => 1, 'ATAC' => 1, 'GCAG' => 1
  }->{$splice};
  return $canonical || 0;
}

=head2 splice_seq

  Example     : my ($donor, $acceptor) = @{$intron->splice_seq}; 
  Description : Get the donor and acceptor splice sites for this intron
  Returntype  : ArrayRef[String] The donor and acceptor sequences as Strings
  Exceptions  : Thrown if a feature Slice cannot be found

=cut

sub splice_seq {
  my ($self) = @_;
  my $slice = $self->feature_Slice();
  throw "Cannot retrieve feature_Slice() for this Intron" unless $slice;
  my $length = $self->length();
  my $donor_seq    = uc($slice->subseq(1,2));
  my $acceptor_seq = uc($slice->subseq($length - 1, $length));
  return [$donor_seq, $acceptor_seq];
}

1;


