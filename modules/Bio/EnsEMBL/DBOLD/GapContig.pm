
#
# Ensembl module for Bio::EnsEMBL::DBOLD::GapContig
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::DBOLD::GapContig - DESCRIPTION of Object

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

Describe the object here

=head1 AUTHOR - Ewan Birney

This modules is part of the Ensembl project http://www.ensembl.org

Email birney@ebi.ac.uk

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...


package Bio::EnsEMBL::DBOLD::GapContig;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::EnsEMBL::Root

use Bio::EnsEMBL::Root;
use Bio::PrimarySeqI;

@ISA = qw(Bio::PrimarySeqI Bio::EnsEMBL::Root);

# new() is written here 

sub new {
  my($class,@args) = @_;

  my $self = {};
  bless $self,$class;
  
# set stuff in self from @args
  return $self;
}


=head2 length

 Title   : length
 Usage   : $obj->length($newval)
 Function: 
 Example : 
 Returns : value of length
 Args    : newvalue (optional)


=cut

sub length{
   my ($obj,$value) = @_;
   return 100000000000;
}

=head2 subseq

 Title   : subseq
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub subseq{
   my ($self,$start,$end) = @_;

   my $len = $end - $start +1;
   return 'N' x $len;
}

=head2 seq

 Title   : seq
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub seq{
   my ($self,@args) = @_;

   $self->throw("Calling seq on GapContig. Bad news!");


}


=head2 id

 Title   : id
 Usage   : $obj->id($newval)
 Function: 
 Example : 
 Returns : value of id
 Args    : newvalue (optional)


=cut

sub id{
   my ($obj,$value) = @_;
   if( defined $value) {
      $obj->{'id'} = $value;
    }
    return $obj->{'id'};

}


=head2 desc

 Title   : desc
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub desc{
   my ($self,@args) = @_;

   return "";
}







