
#
# BioPerl module for Bio::EnsEMBL::Homol
#
# Cared for by Ewan Birney <birney@sanger.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Homol - DESCRIPTION of Object

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


package Bio::EnsEMBL::Homol;
use vars qw(@ISA);
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::SeqFeature::Homol;

@ISA = qw(Bio::SeqFeature::Homol);
# new() is inherited from Bio::Root::Object

# _initialize is where the heavy stuff will happen when new is called

sub _initialize {
  my($self,@args) = @_;

  my $make = $self->SUPER::_initialize;

# set stuff in self from @args
 return $make; # success - we hope!
}

=head2 analysis

 Title   : analysis
 Usage   : $sf->analysis();
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub analysis{
   my ($self,$value) = @_;

   # This is just so we are in sync with Michele's stuff.
   # should not occur in the future! 
   # FIXME: MUST REMOVE ASAP
   if( defined $value ) {
       if( ! ref $value || !$value->isa('Bio::EnsEMBL::Analysis::Analysis') ) {
	   $self->throw("Trying to add a non analysis object! [$value]");
       }
       # this should be stored in a different way
       $self->add_tag_value('Analysis',$value);
   }

   my @ref = $self->each_tag_value('Analysis');
   return $ref[0];
   
}


=head2 to_FTHelper

 Title   : to_FTHelper
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub to_FTHelper{
   my ($self) = @_;

   # Make new FTHelper, and fill in the key
   my $fth = Bio::AnnSeqIO::FTHelper->new;
   $fth->key('similarity');
   
   # Add location line
   my $g_start = $self->start;
   my $g_end   = $self->end;
   my $loc = "$g_start..$g_end";
   if ($self->strand == -1) {
        $loc = "complement($loc)";
    }
   $fth->loc($loc);
   
   # Add note describing similarity
   my $type    = $self->homol_SeqFeature->seqname;
   my $r_start = $self->homol_SeqFeature->start;
   my $r_end   = $self->homol_SeqFeature->end;
   $fth->add_field('note', "$type: matches $r_start to $r_end");
   $fth->add_field('note', "score=".$self->score);
   
   
   return $fth;
}
1;
