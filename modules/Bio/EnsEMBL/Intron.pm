#
# EnsEMBL module for Bio::EnsEMBL::Exon
#
#
# Copyright Ian Longden
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=pod 

=head1 NAME Bio::EnsEMBL::Intron - A class representing an Intron

=head1 SYNOPSIS

    $ex = new Bio::EnsEMBL::Intron(exon e1, exon e2);

=cut


package Bio::EnsEMBL::Intron;
use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::SeqFeature;
use Bio::Seq; # introns have to have sequences...

@ISA = qw(Bio::EnsEMBL::Root Bio::EnsEMBL::SeqFeature);

=head2 new

  Args       : exon1, exon2. The two exons to build the Intron from.
  Example    : $intron = new Bio::EnsEMBL::Intron($exon1, $exon2)
  Description: create an Intron object
  Returntype : Bio::EnsEMBL::Intron
  Exceptions : exons not on the same strand or slice.
  Caller     : general

=cut

sub new {
  my ($class,$e1,$e2) = @_;

  $class = ref $class || $class;

  my $self = $class->SUPER::new();

  if($e1->strand == -1){
    $self->end(($e1->start)-1);
    $self->start(($e2->end)+1);
  }
  else{
    $self->start(($e1->end)+1);
    $self->end(($e2->start)-1);
  }

  if($e1->strand != $e2->strand){
    $self->throw("Exons on different strand. Not allowed");
  }
  else{
    $self->strand($e1->strand);
  }

  if($e1->contig != $e2->contig){
    $self->throw("Exons on different slices. Not allowed");
  }
  else{ 
    $self->contig($e1->contig);
  }

  $self->prev_Exon($e1);
  $self->next_Exon($e2);

  return $self;
}


=head2 prev_Exon

  Args       : none
  Example    : $exon = $intron->prev_Exon
  Description: Returns the exon before this Intron
  Returntype : Bio::EnsEMBL::Exon
  Exceptions : none
  Caller     : general

=cut


sub prev_Exon {
    my( $self, $prev_Exon ) = @_;
    
    if ($prev_Exon) {
        $self->{'_prev_Exon'} = $prev_Exon;
    }
    return $self->{'_prev_Exon'};
}


=head2 next_Exon

  Args       : none
  Example    : $exon = $intron->next_Exon
  Description: Returns the exon after this Intron
  Returntype : Bio::EnsEMBL::Exon
  Exceptions : none
  Caller     : general

=cut

sub next_Exon {
    my( $self, $next_Exon ) = @_;
    
    if ($next_Exon) {
        $self->{'_next_Exon'} = $next_Exon;
    }
    return $self->{'_next_Exon'};
}

sub stable_id {
    my( $self ) = @_;
    
    my $prev = $self->prev_Exon;
    my $next = $self->next_Exon;
    return $self->prev_Exon->stable_id
        . "-"
        . $self->next_Exon->stable_id;
}

1;


