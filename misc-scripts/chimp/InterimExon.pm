use strict;
use warnings;

use StatMsg;

package InterimExon;

#
# errors which are fatal for exons
#
my @FATAL =
  (StatMsg::DELETE | StatMsg::CDS | StatMsg::LONG,
   StatMsg::DELETE | StatMsg::CDS | StatMsg::MEDIUM | StatMsg::FRAMESHIFT,
   StatMsg::STRAND_FLIP,
   StatMsg::INVERT,
   StatMsg::SCAFFOLD_SPAN,
   StatMsg::CONFUSED);


sub new {
  my $class = shift;

  return bless {'StatMsgs' => [],
                'fail'     => 0}, $class;
}


#
# returns true if this exon has a 'fatal' error
#
sub is_fatal {
  my $self = shift;

  foreach my $msg (@{$self->get_all_StatMsgs}) {
    foreach my $code (@FATAL) {
      if(($msg & $code) == $code) {
	return 1;
      }
    }
  }

  return 0;
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


sub start_phase {
  my $self = shift;
  $self->{'start_phase'} = shift if(@_);
  return $self->{'start_phase'};
}


sub end_phase {
  my $self = shift;
  $self->{'end_phase'} = shift if(@_);
  return $self->{'end_phase'};
}

sub fail {
  my $self = shift;
  $self->{'fail'} = shift if(@_);
  return $self->{'fail'};
}





1;
