
#
# BioPerl module for ContigOverlap
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::ContigOverlap - Describes the overlap between two contigs

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

The contig overlap object is directional - you pick it up from a particular
contig and it tells you what position on that contig is equivalent to what
position on the other contig.

http://www.ncbi.nlm.nih.gov/genome/seq/Hs_Data/contig.xml

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::ContigOverlap;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::Object;


@ISA = qw(Bio::Root::Object);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize(@args);

  # set stuff in self from @args
  return $make; # success - we hope!
}

=head2 sister_contig_id

 Title   : sister_contig_id
 Usage   : $obj->sister_contig_id($newval)
 Function: 
 Returns : value of sister_contig_id
 Args    : newvalue (optional)


=cut

sub sister_contig_id{
   my $obj = shift;
   if( @_ ) {
       my $value = shift;
       $obj->{'sister_contig_id'} = $value;
   }
   return $obj->{'sister_contig_id'};
   
}

=head2 sister_contig_position
    
 Title   : sister_contig_position
 Usage   : $obj->sister_contig_position($newval)
 Function: 
 Returns : value of sister_contig_position
 Args    : newvalue (optional)


=cut

sub sister_contig_position{
   my $obj = shift;
   if( @_ ) {
       my $value = shift;
       $obj->{'sister_contig_position'} = $value;
   }
   return $obj->{'sister_contig_position'};
   
}

=head2 sister_contig_polarity

 Title   : sister_contig_polarity
 Usage   : $obj->sister_contig_polarity($newval)
 Function: 
 Returns : value of sister_contig_polarity
 Args    : newvalue (optional)


=cut

sub sister_contig_polarity {
   my $obj = shift;
   if( @_ ) {
       my $value = shift;
       $obj->{'sister_contig_polarity'} = $value;
   }
   return $obj->{'sister_contig_polarity'};
   
}

=head2 self_contig_position

 Title   : self_contig_position
 Usage   : $obj->self_contig_position($newval)
 Function: 
 Returns : value of self_contig_position
 Args    : newvalue (optional)


=cut

sub self_contig_position{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'self_contig_position'} = $value;
  }
   return $obj->{'self_contig_position'};
   
}


1;
