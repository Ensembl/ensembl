
#
# BioPerl module for Bio::EnsEMBL::AnnSeq
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::AnnSeq - EnsEMBL Annotated Sequence - derived from Bioperl AnnSeq

=head1 SYNOPSIS

Give standard usage here

=head1 DESCRIPTION

EnsEMBL requires a number of specialised fields in the AnnSeq object so that
we can dump out the correct information etc. Rather than making a pig's ear
of the Bioperl AnnSeq object, we are factoring in this stuff ontop here.

This object is-a AnnSeq object, and the core functionality is there. Read
the Bio::AnnSeq notes first.

=head1 CONTACT

Describe contact details here

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut


# Let the code begin...'


package Bio::EnsEMBL::AnnSeq;
use vars qw($AUTOLOAD @ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Seq;

@ISA = qw(Bio::Seq);

# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize(@args);

# set stuff in self from @args
  return $make; # success - we hope!
}

=head2 htg_phase

 Title   : htg_phase
 Usage   : $obj->htg_phase($newval)
 Function: 
 Returns : value of htg_phase
 Args    : newvalue (optional)


=cut

sub htg_phase{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'htg_phase'} = $value;
    }
    return $obj->{'htg_phase'};

}


=head2 sv

 Title   : sv
 Usage   : $obj->sv($newval)
 Function: 
 Returns : value of sv
 Args    : newvalue (optional)


=cut

sub sv{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'sv'} = $value;
    }
    return $obj->{'sv'};

}


=head2 embl_id

 Title   : embl_id
 Usage   : $obj->embl_id($newval)
 Function: 
 Returns : value of embl_id
 Args    : newvalue (optional)


=cut

sub embl_id{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'embl_id'} = $value;
    }
    return $obj->{'embl_id'};

}

=head2 project_name

 Title   : project_name
 Usage   : $obj->project_name($newvalue)
 Function: Stores the name of the Sanger Centre project.
           Used only for generating the special "AC * "
           lines for dumping EMBL files suitable for
           submission to EMBL.
           See: Bio::EnsEMBL::EMBL_Dump_Sanger
 Returns : value of project_name
 Args    : newvalue (optional)


=cut

sub project_name {
    my $obj = shift;
    if( @_ ) {
        $obj->{'project_name'} = shift;
    }
    return $obj->{'project_name'};
}
