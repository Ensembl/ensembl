
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

use POSIX;
use strict;


use Bio::EnsEMBL::Utils::Exception qw(warning);
BEGIN {
 eval {
 require Time::HiRes;
 Time::HiRes->import('time');
 };
};


sub new { 
    my ($class,$name) = @_;
    my $self = {
       'is_active'       => 0,
       'total_time'      => 0,
       'total_time_time' => 0,
       'max_time'        => 0, 
       'min_time'        => 999999999,
       'number'          => 0,
       'tag'             => $name
    };
    bless $self,$class;
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

   if( $self->{'is_active'} == 1 ) {
       warning("Attempting to push stack on tag ".$self->tag." when active. Discarding previous push");
   }
   #my($user,$sys) = times();
   # $self->{'current_start'} = (POSIX::times)[0];
   $self->{'current_start'} = time();
   $self->{'is_active'}=1
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

   if( $self->{'is_active'} == 0 ) {
       warning("Attempting to pop stack on tag ".$self->tag." when not active. Ignoring");
   }
   #my($user,$sys) = times();
 #  my $clocktime = ( (POSIX::times)[0] - $self->{'current_start'} ) / POSIX::sysconf(&POSIX::_SC_CLK_TCK);
   my $clocktime = time() - $self->{'current_start'};
   $self->{'max_time'} = $clocktime if $self->{'max_time'} < $clocktime;
   $self->{'min_time'} = $clocktime if $self->{'min_time'} > $clocktime;
   $self->{'total_time'}+=$clocktime;
   $self->{'total_time_time'} += $clocktime * $clocktime;
   $self->{'number'}++;
   $self->{'is_active'}=0;
}


=head2 total_time_time

 Title   : total_time_time
 Usage   : $obj->total_time_time($newval)
 Function: 
 Returns : value of total_time_time
 Args    : newvalue (optional)


=cut

sub total_time_time {
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'total_time_time'} = $value;
    }
    return $obj->{'total_time_time'};

}

=head2 max_time

 Title   : max_time
 Usage   : $obj->max_time($newval)
 Function: 
 Returns : value of max_time
 Args    : newvalue (optional)


=cut

sub max_time{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'max_time'} = $value;
    }
    return $obj->{'max_time'};
}

=head2 min_time

 Title   : min_time
 Usage   : $obj->min_time($newval)
 Function: 
 Returns : value of min_time
 Args    : newvalue (optional)


=cut

sub min_time{
   my $obj = shift;
   if( @_ ) {
      my $value = shift;
      $obj->{'min_time'} = $value;
    }
    return $obj->{'min_time'};
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
