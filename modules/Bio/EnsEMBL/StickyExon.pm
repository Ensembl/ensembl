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

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

  # Array to store exon tags
  $self->{_component_exons} = [];

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

1;
