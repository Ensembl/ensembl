use strict;
use warnings;

package InterimTranscript;

use Bio::EnsEMBL::Utils::Exception qw(warning);



sub new {
  my $class = shift;

  return bless {'exons' => [],
                'StatsMsgs' => []}, $class;
}


sub add_StatMsg {
  my $self = shift;
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


sub add_Exon {
  my $self = shift;
  my $exon = shift;

  push @{$self->{'exons'}}, $exon;
}

sub get_all_Exons {
  my $self = shift;

  return $self->{'exons'};
}

sub flush_Exons {
  my $self = shift;
  $self->{'exons'} = [];
}


sub stable_id {
  my $self = shift;
  $self->{'stable_id'} = shift if(@_);
  return $self->{'stable_id'};
}

sub version {
  my $self = shift;
  $self->{'version'} = shift if(@_);
  return $self->{'version'};
}

sub cdna_coding_start {
  my $self = shift;
  $self->{'cdna_coding_start'} = shift if(@_);
  return $self->{'cdna_coding_start'};
}

sub cdna_coding_end {
  my $self = shift;
  $self->{'cdna_coding_end'} = shift if(@_);
  return $self->{'cdna_coding_end'};
}


sub move_cdna_coding_start {
  my $self = shift;
  my $offset = shift;
  $self->{'cdna_coding_start'} += $offset;
}

sub move_cdna_coding_end {
  my $self = shift;
  my $offset = shift;
  $self->{'cdna_coding_end'} += $offset;
}


1;
