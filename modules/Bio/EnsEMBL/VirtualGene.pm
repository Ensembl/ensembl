
#
# BioPerl module for Bio::EnsEMBL::VirtualGene
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::VirtualGene - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::VirtualGene;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::SeqFeatureI;
use Bio::Root::Object

@ISA = qw(Bio::Root::Object Bio::SeqFeatureI);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize(@args);

  my ($gene,$contig) = $self->_rearrange(['GENE','CONTIG'],@args);
  if( !defined $gene ) {
      $self->throw("No gene in virtualgene object");
  }
  if( !defined $contig || ! ref $contig || ! $contig->isa('Bio::EnsEMBL::DB::ContigI') ) {
      $self->throw("you have to have a virtual gene on a particular contig");
  }


  $self->gene($gene);

  return $make; # success - we hope!
}

=head2 start

 Title   : start
 Usage   : $obj->start($newval)
 Function: 
 Returns : value of start
 Args    : newvalue (optional)


=cut

sub start{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'start'} = $value;
    }
    return $obj->{'start'};

}

=head2 end

 Title   : end
 Usage   : $obj->end($newval)
 Function: 
 Returns : value of end
 Args    : newvalue (optional)


=cut

sub end{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'end'} = $value;
    }
    return $obj->{'end'};

}

=head2 strand

 Title   : strand
 Usage   : $obj->strand($newval)
 Function: 
 Returns : value of strand
 Args    : newvalue (optional)


=cut

sub strand{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'strand'} = $value;
    }
    return $obj->{'strand'};

}

=head2 seqname

 Title   : seqname
 Usage   : $obj->seqname($newval)
 Function: 
 Returns : value of seqname
 Args    : newvalue (optional)


=cut

sub seqname{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'seqname'} = $value;
    }
    return $obj->{'seqname'};

}



=head2 gene

 Title   : gene
 Usage   : $obj->gene($newval)
 Function: 
 Returns : value of gene
 Args    : newvalue (optional)


=cut

sub gene{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      if( ! ref $value || ! $value->isa("Bio::EnsEMBL::Gene") ) {
	  $self->throw("Gene object must inheriet from Gene...");
      }

      $obj->{'gene'} = $value;
    }
    return $obj->{'gene'};

}

1;
