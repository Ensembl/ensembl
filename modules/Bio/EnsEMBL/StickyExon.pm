#
# BioPerl module for Exon
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=pod 

=head1 NAME

Bio::EnsEMBL::StickyExon - A Confirmed Exon which spans two contigs internally

=head1 SYNOPSIS

    $sticky = new Bio::EnsEMBL::Exon;

    # is a normal exon
    $sticky->start();
    $sticky->end();

    # has component_Exons
    foreach $sub ( $sticky->each_component_Exon ) {
       # $sub is an exon that ends on a contig
    }

=head1 DESCRIPTION

Sticky Exons represent Exons which internally span contigs. They are made during the
write back on virtual contigs, which writes the exons that span joins into the database.


=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::StickyExon;
use vars qw(@ISA $AUTOLOAD);
use strict;

# Object preamble - inherits from Bio::SeqFeature::Generic

use Bio::EnsEMBL::Exon;


@ISA = qw(Bio::EnsEMBL::Exon);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub new {
  my($class,@args) = @_;
  
  my $self = Bio::EnsEMBL::Exon->new(@args);
  bless $self,$class;



  # Array to store exon tags
  $self->{_component_exons} = [];
  
  return $self;
}


=head2 id

 Title   : id
 Usage   : overrides id to get/set locally and sets all component exons
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub id{
   my ($self,$value) = @_;

   if( defined $value ) {
       $self->{'_sticky_id'} = $value;
       foreach my $c ( $self->each_component_Exon() ) {
	   $c->id($value);
       }
   }

   return $self->{'_sticky_id'};

}


=head2 each_component_Exon

 Title   : each_component_Exon
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub each_component_Exon{
   my ($self,@args) = @_;

   return @{$self->{'_component_exons'}};
}



=head1

  Arg  1   : integer start - relative to the exon
  Arg  2   : integer end   - relative to the exon

  Function : Provides a list of Bio::EnsEMBL::SeqFeatures which
             is the genomic coordinates of this start/end on the exon
             For simple exons this is one feature  for Stickies this
             is overridden and gives out a list of Bio::EnsEMBL::SeqFeatures

  Returns  : list of Bio::EnsEMBL::SeqFeature


=cut

sub contig_seqfeatures_from_relative_position {
  my ($self,$start,$end) = @_;

  if( !defined $end ) {
    $self->throw("Have not defined all the methods!");
  }

  # easy
  if( $start < 1 ) {
    $self->warn("Attempting to fetch start less than 1 ($start)");
    $start = 1;
  }

  if( $end > $self->length ) {
    $self->warn("Attempting to fetch end greater than end of exon ($end)");
    $end = $self->length;
  }

  my @out;
  my $sf;
  my @exons = $self->each_component_Exon();
  my $len = 0;
  while( scalar(@exons) > 0 ) {
    if( $exons[0]->length + $len > $start ) {
       last;
    } else {
       my $discard = shift @exons;
       $len += $discard;
    }
  }

  # handle the first component exon

  if( scalar(@exons) == 0 ) {
     return @out;
  }
  
  $sf = Bio::EnsEMBL::SeqFeature->new();
  $sf->seqname($exons[0]->contig->id);
  $sf->strand($exons[0]->strand);
  $sf->start($exons[0]->start + $start - $len);

  if( $end < $len + $exons[0]->length ) {
      $sf->end($exons[0]->start + $end - $len);
      return $sf;
  } else {
      $sf->end($exons[0]->end);
      push(@out,$sf);
  }


  while( scalar(@exons) ) {
     if( $exons[0]->length + $len > $end ) {
        last;
     }
     $sf = Bio::EnsEMBL::SeqFeature->new();
     $sf->seqname($exons[0]->contig->id);
     $sf->strand($exons[0]->strand);
     $sf->start($exons[0]->start);
     $sf->start($exons[0]->end);
     push(@out,$sf);
     $len += $exons[0]->length;
  }

  if( scalar(@exons) == 0 ) {
     return @out;
  }

  # handle the last exon

  $sf = Bio::EnsEMBL::SeqFeature->new();
  $sf->seqname($exons[0]->contig->id);
  $sf->strand($exons[0]->strand);
  $sf->start($exons[0]->start);
  $sf->start($exons[0]->start + $end - $len);



  return @out;
}

=head2 add_component_Exon

 Title   : add_component_Exon
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_component_Exon{
   my ($self,$exon) = @_;

   if( !ref $exon || ! $exon->isa('Bio::EnsEMBL::Exon')) {
       $self->throw("$exon is not an exon");
   }

   push(@{$self->{'_component_exons'}},$exon);
}

=head2 length

 Title   : length
 Usage   : length
 Function: calculate number of  nucleotides constituting the Exon
 Example :
 Returns : a number
 Args    : none

=cut


sub length {
    my $self = shift;

    my $len =0; 

    foreach my $subexon ( $self->each_component_Exon ) {
        $len += $subexon->length;
    }
    return $len;
}

=head2 _sort_by_sticky_rank

 Title   : _sort_by_sticky_rank
 Usage   : 
 Function: put the contained exons in the right order
 Example :
 Returns : 
 Args    : 

=cut


sub _sort_by_sticky_rank {
    my $self = shift;

    my @sorted;

    @sorted= sort {$a->sticky_rank <=> $b->sticky_rank } 
      @{$self->{'_component_exons'}};
    $self->{'_component_exons'} = \@sorted;
    return 1;
}



1;
