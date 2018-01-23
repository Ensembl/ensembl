=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::Util::EprofStack - DESCRIPTION of Object

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 METHODS

=cut

package Bio::EnsEMBL::Utils::EprofStack;

use strict;
use warnings;

use POSIX;

use Bio::EnsEMBL::Utils::Exception ('warning');

BEGIN {
  eval {
    require Time::HiRes;
    Time::HiRes->import('time');
  };
}

sub new {
  my ( $proto, $name ) = @_;

  my $class = ref($proto) || $proto;

  my $self = bless( { 'is_active'       => 0,
                      'total_time'      => 0,
                      'total_time_time' => 0,
                      'max_time'        => 0,
                      'min_time'        => 999999999,
                      'number'          => 0,
                      'tag'             => $name
                    },
                    $class );

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

sub push_stack {
  my ( $self, @args ) = @_;

  if ( $self->{'is_active'} == 1 ) {
    warning(
             sprintf(     "Attempting to push stack on tag '%s' "
                        . "when active. Discarding previous push."
                        . $self->tag() ) );
  }

  # my ( $user, $sys ) = times();
  # $self->{'current_start'} = (POSIX::times)[0];

  $self->{'current_start'} = time();
  $self->{'is_active'}     = 1;
}

=head2 pop_stack

 Title   : pop_stack
 Usage   :
 Function:
 Example :
 Returns : 
 Args    :


=cut

sub pop_stack {
  my ( $self, @args ) = @_;

  if ( $self->{'is_active'} == 0 ) {
    warning(
             sprintf( "Attempting to pop stack on tag '%s' "
                        . "when not active. Ignoring.",
                      $self->tag() ) );
  }

  # my ( $user, $sys ) = times();
  # my $clocktime =
  #   ( (POSIX::times)[0] - $self->{'current_start'} )/
  #   POSIX::sysconf(&POSIX::_SC_CLK_TCK);

  my $clocktime = time() - $self->{'current_start'};

  if ( $self->{'max_time'} < $clocktime ) {
    $self->{'max_time'} = $clocktime;
  }
  if ( $self->{'min_time'} > $clocktime ) {
    $self->{'min_time'} = $clocktime;
  }

  $self->{'total_time'}      += $clocktime;
  $self->{'total_time_time'} += $clocktime*$clocktime;
  $self->{'number'}++;
  $self->{'is_active'} = 0;
} ## end sub pop_stack

=head2 total_time_time

 Title   : total_time_time
 Usage   : $obj->total_time_time($newval)
 Function: 
 Returns : value of total_time_time
 Args    : newvalue (optional)


=cut

sub total_time_time {
  my ( $self, $value ) = @_;

  if ( defined($value) ) { $self->{'total_time_time'} = $value }

  return $self->{'total_time_time'};
}

=head2 max_time

 Title   : max_time
 Usage   : $obj->max_time($newval)
 Function: 
 Returns : value of max_time
 Args    : newvalue (optional)


=cut

sub max_time {
  my ( $self, $value ) = @_;

  if ( defined($value) ) { $self->{'max_time'} = $value }

  return $self->{'max_time'};
}

=head2 min_time

 Title   : min_time
 Usage   : $obj->min_time($newval)
 Function: 
 Returns : value of min_time
 Args    : newvalue (optional)


=cut

sub min_time {
  my ( $self, $value ) = @_;

  if ( defined($value) ) { $self->{'min_time'} = $value }

  return $self->{'min_time'};
}

=head2 total_time

 Title   : total_time
 Usage   : $obj->total_time($newval)
 Function: 
 Returns : value of total_time
 Args    : newvalue (optional)


=cut

sub total_time {
  my ( $self, $value ) = @_;

  if ( defined($value) ) { $self->{'total_time'} = $value }

  return $self->{'total_time'};
}

=head2 number

 Title   : number
 Usage   : $obj->number($newval)
 Function: 
 Returns : value of number
 Args    : newvalue (optional)


=cut

sub number {
  my ( $self, $value ) = @_;

  if ( defined($value) ) { $self->{'number'} = $value }

  return $self->{'number'};
}

=head2 is_active

 Title   : is_active
 Usage   : $obj->is_active($newval)
 Function: 
 Returns : value of is_active
 Args    : newvalue (optional)


=cut

sub is_active {
  my ( $self, $value ) = @_;

  if ( defined($value) ) { $self->{'is_active'} = $value }

  return $self->{'is_active'};
}

=head2 current_start

 Title   : current_start
 Usage   : $obj->current_start($newval)
 Function: 
 Returns : value of current_start
 Args    : newvalue (optional)


=cut

sub current_start {
  my ( $self, $value ) = @_;

  if ( defined($value) ) { $self->{'current_start'} = $value }

  return $self->{'current_start'};
}

=head2 tag

 Title   : tag
 Usage   : $obj->tag($newval)
 Function: 
 Returns : value of tag
 Args    : newvalue (optional)


=cut

sub tag {
  my ( $self, $value ) = @_;

  if ( defined($value) ) { $self->{'tag'} = $value }

  return $self->{'tag'};
}

1;
