
#
# BioPerl module for Bio::EnsEMBL::Util::EprofStack
#
# Cared for by Ewan Birney <birney@ebi.ac.uk>
#
# Copyright Ewan Birney
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Util::EprofStack - DESCRIPTION of Object

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


package Bio::EnsEMBL::Utils::EprofStack;
use vars qw(@ISA);
use POSIX;
use strict;

# Object preamble - inheriets from Bio::Root::Object

use Bio::Root::RootI;

@ISA = qw(Bio::Root::RootI);


sub new { 
    my ($class,$name) = @_;
    my $self = {};
    bless $self,$class;

    $self->is_active(0);
    $self->total_time(0);
    $self->number(0);
    $self->tag($name);

    return $self;
}

=head2 push_stack

 Title   : push_stack
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub push_stack{
   my ($self,@args) = @_;

   if( $self->is_active == 1 ) {
       $self->warn("Attempting to push stack on tag ",$self->tag," when active. Discarding previous push");
   }
   #my($user,$sys) = times();
   my $real = (POSIX::times)[0];
   $self->current_start($real);
   $self->is_active(1);
}

=head2 pop_stack

 Title   : pop_stack
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub pop_stack{
   my ($self,@args) = @_;

   if( $self->is_active == 0 ) {
       $self->warn("Attempting to push stack on tag ",$self->tag," when not active. Ignoring");
   }
   #my($user,$sys) = times();
   my $real = (POSIX::times)[0];
   my $clockticks =  $real - $self->current_start;
   my $clocktime = $clockticks / POSIX::sysconf(&POSIX::_SC_CLK_TCK);
   $self->total_time( $self->total_time() + $clocktime);
   $self->number( $self->number() + 1 );
   $self->is_active(0);
}


=head2 total_time

 Title   : total_time
 Usage   : $obj->total_time($newval)
 Function: 
 Returns : value of total_time
 Args    : newvalue (optional)


=cut

sub total_time{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'total_time'} = $value;
    }
    return $obj->{'total_time'};

}

=head2 number

 Title   : number
 Usage   : $obj->number($newval)
 Function: 
 Returns : value of number
 Args    : newvalue (optional)


=cut

sub number{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'number'} = $value;
    }
    return $obj->{'number'};

}

=head2 is_active

 Title   : is_active
 Usage   : $obj->is_active($newval)
 Function: 
 Returns : value of is_active
 Args    : newvalue (optional)


=cut

sub is_active{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'is_active'} = $value;
    }
    return $obj->{'is_active'};

}

=head2 current_start

 Title   : current_start
 Usage   : $obj->current_start($newval)
 Function: 
 Returns : value of current_start
 Args    : newvalue (optional)


=cut

sub current_start{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'current_start'} = $value;
    }
    return $obj->{'current_start'};

}


=head2 tag

 Title   : tag
 Usage   : $obj->tag($newval)
 Function: 
 Returns : value of tag
 Args    : newvalue (optional)


=cut

sub tag{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'tag'} = $value;
    }
    return $obj->{'tag'};

}

1;
