
#
# Ensembl module for Bio::EnsEMBL::Mapper::Pair
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Mapper::Pair - DESCRIPTION of Object

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


package Bio::EnsEMBL::Mapper::Pair;
use vars qw(@ISA);
use strict;

# Object preamble - inherits from Bio::Root::RootI

use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI);

# new() is written here 

sub new {
  my($class,@args) = @_;

    my $self = {};
    bless $self,$class;

# set stuff in self from @args
    return $self;
}

=head2 to

 Title   : to
 Usage   : $obj->to($newval)
 Function: 
 Example : 
 Returns : value of to
 Args    : newvalue (optional)


=cut

sub to{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'to'} = $value;
    }
    return $self->{'to'};

}

=head2 from

 Title   : from
 Usage   : $obj->from($newval)
 Function: 
 Example : 
 Returns : value of from
 Args    : newvalue (optional)


=cut

sub from{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'from'} = $value;
    }
    return $self->{'from'};

}

=head2 ori

 Title   : ori
 Usage   : $obj->ori($newval)
 Function: 
 Example : 
 Returns : value of ori
 Args    : newvalue (optional)


=cut

sub ori{
   my ($self,$value) = @_;
   if( defined $value) {
      $self->{'ori'} = $value;
    }
    return $self->{'ori'};

}

1;
