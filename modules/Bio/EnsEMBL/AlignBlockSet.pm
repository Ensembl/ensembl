
#
# Ensembl module for Bio::EnsEMBL::AlignBlockSet
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::AlignBlockSet - DESCRIPTION of Object

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


package Bio::EnsEMBL::AlignBlockSet;
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
    $self->{'_align_block'}= [];
# set stuff in self from @args
    return $self;
}

=head2 get_AlignBlocks

 Title   : get_AlignBlocks
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub get_AlignBlocks{
   my ($self,@args) = @_;

   return @{$self->{'_align_block'}};
}


=head2 add_AlignBlock

 Title   : add_AlignBlock
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub add_AlignBlock{
   my ($self,$bl) = @_;

   if( !ref $bl || !$bl->isa('Bio::EnsEMBL::AlignBlock') ) {
       $self->throw("Must add an AlignBlock [$bl]");
   }

   push(@{$self->{'_align_block'}},$bl);
}

1;



