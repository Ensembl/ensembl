use strict;
use warnings;

package InterimExon;


sub new {
  my $class = shift;

  return bless {'StatMsgs' => [],
                'fail'     => 0}, $class;
}

sub add_StatMsg {
  my $self    = shift;
  my $statMsg = shift;
  push @{$self->{'StatMsgs'}}, $statMsg;
}

sub get_all_StatMsgs {
  my $self = shift;
  return @{$self->{'StatMsgs'}};
}

sub last_StatMsg {
  my $self = shift;

  my @msgs = @{$self->{'StatMsgs'}};
  return undef if(!@msgs);
  return $msgs[$#msgs];
}

sub flush_StatMsgs {
  my $self = shift;
  $self->{'StatMsgs'} = [];
}

sub stable_id {
  my $self = shift;
  $self->{'stable_id'} = shift if(@_);
  return $self->{'stable_id'};
}

sub start {
  my $self = shift;
  $self->{'start'} = shift if(@_);
  return $self->{'start'};
}

sub end {
  my $self = shift;
  $self->{'end'} = shift if(@_);
  return $self->{'end'};
}

sub length {
  my $self = shift;
  return $self->end() - $self->start() + 1;
}

sub strand {
  my $self = shift;
  $self->{'strand'} = shift if(@_);
  return $self->{'strand'};
}

sub seq_region {
  my $self = shift;
  $self->{'seq_region'} = shift if(@_);
  return $self->{'seq_region'};
}

sub cdna_start {
  my $self = shift;
  $self->{'cdna_start'} = shift if(@_);
  return $self->{'cdna_start'};
}


sub cdna_end {
  my $self = shift;
  $self->{'cdna_end'} = shift if(@_);
  return $self->{'cdna_end'};
}


sub slice {
  my $self = shift;
  $self->{'slice'} = shift if(@_);
  return $self->{'slice'};
}

sub fail {
  my $self = shift;
  $self->{'fail'} = shift if(@_);
  return $self->{'fail'};
}





1;
