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

=head2 add_Supporting_Feature

 Title   : add_Supporting_Feature
 Usage   : $obj->add_Supporting_Feature($feature)
 Function: 
 Returns : Nothing
 Args    : Bio::EnsEMBL::SeqFeature


=cut


sub add_Supporting_Feature {
    my ($self,$feature) = @_;

    $self->throw("Supporting evidence [$feature] not Bio::EnsEMBL::SeqFeatureI") unless 
	defined($feature) &&  $feature->isa("Bio::EnsEMBL::SeqFeatureI");

    $self->{_supporting_evidence} = undef unless defined($self->{_supporting_evidence});

    # check whether this feature object has been added already
    my $found = 0;
    if ( $feature && $self->{_supporting_evidence} ){
      foreach my $added_feature ( @{ $self->{_supporting_evidence} } ){
	# compare objects
	if ( $feature == $added_feature ){
	  $found = 1;
	  
	  # no need to look further
	  last;
	}
      }
    }
    if ( $found == 0 ){
      push(@{$self->{_supporting_evidence}},$feature);
    }
}



=head2 each_Supporting_Feature

 Title   : each_Supporting_Feature
 Usage   : my @f = $obj->each_Supporting_Feature
 Function: 
 Returns : @Bio::EnsEMBL::Feature
 Args    : none


=cut


sub each_Supporting_Feature {
    my ($self) = @_;

    if ( !defined ( $self->{_supporting_evidence} )) {
      $self->{_supporting_evidence} = [];  
      $self->adaptor->fetch_evidence_by_Exon( $self );
    }

    return @{$self->{_supporting_evidence}};
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

=head1 cdna_coord_2_features

  Arg  1   : integer start - relative to the exon
  Arg  2   : integer end   - relative to the exon

  Function : Provides a list of Bio::EnsEMBL::SeqFeatures which
             is the genomic coordinates of this start/end on the exon
             For simple exons this is one feature - for Stickies this
             is overridden

  Returns  : list of Bio::EnsEMBL::SeqFeature


=cut

sub cdna_coord_2_features {
  my ($self,$start,$end) = @_;

  my ( @features, @result );
  # ec start end are cdna relative positions for the current exons cdna
  my ( $ec_start, $ec_end, $offset );

  # which area of the exon overlaps with requested start,end
  # usually either the last bit, the first bit or all.
  my ( $ov_start, $ov_end );

  $ec_start = 0; $ec_end = 0;

  my @exons = $self->each_component_Exon();

  while( my $exon = shift  @exons ) {

    if( $ec_start > 0 ) {
      $ec_start = $ec_end + 1;
      $ec_end = $exon->length() + $ec_start - 1;
    } else {
      $ec_start = 1;
      $ec_end = $exon->length()
    }
    
    # now exon covers ec_start - ec_end in cdna
    $ov_start = ( $start >= $ec_start ) ? $start : $ec_start;
    $ov_end = ( $end <= $ec_end ) ? $end : $ec_end;
    
    if( $ov_end >= $ov_start ) {
      @features = $exon->cdna_coord_2_features
	( $ov_start-$ec_start + 1, 
	  $ov_end-$ec_start + 1 );
      push( @result, @features );
    }
  }
  return @result;
}

1;
