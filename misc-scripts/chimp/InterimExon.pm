use strict;
use warnings;


package InterimExon;

use Bio::EnsEMBL::Utils::Exception qw(info warning throw);

use StatMsg;

#
# errors which are fatal for exons
#
my @FATAL =
  (StatMsg::DELETE | StatMsg::CDS | StatMsg::LONG,
   StatMsg::INSERT | StatMsg::CDS | StatMsg::LONG,
   StatMsg::DELETE | StatMsg::CDS | StatMsg::MEDIUM | StatMsg::FRAMESHIFT,
   StatMsg::INSERT | StatMsg::CDS | StatMsg::MEDIUM | StatMsg::FRAMESHIFT,
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
      if(($msg->code() & $code) == $code) {
        #info("Code is Fatal: ". StatMsg::code2str($msg->code()));
        return 1;
      }
    }
    #info("Code is NON fatal=". StatMsg::code2str($msg->code()));
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
  return $self->{'StatMsgs'};
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

  if(@_) {
    my $fail = shift;
    #warning("Setting ".$self->stable_id." to failed.\n") if($fail);
    $self->{'fail'} = $fail;
  }

  return $self->{'fail'};
}


#
# This fixes the start and end phases which can get messed up when the
# UTR and CDS move a bit.  We maintain correct phase throughout the program
# but there are a couple of problems that may arise:
#   * UTR at one end of an exon may be completely deleted so the start/end phase
#     needs to change from -1 to 0.
#   * CDS may shrink due to a deletion at 5prime or 3prime end of CDS. This
#     that exons which had start or end phase to have a start phase of -1 now.

sub fix_phase {
  my $exon = shift;
  my $transcript = shift;

  # do not deal with failed exons, we have no chimp coords for these anyway
  # since they were completely lost.
  return if($exon->fail());

  if(!defined($transcript->cdna_coding_start)) {
    throw("cdna coding start not defined!");
  }
  if(!defined($transcript->cdna_coding_end)) {
    throw("cdna coding end not defined.");
  }

  if(!defined($exon->cdna_end())) {
    throw("exons cdna coding end not defined.");
  }

  if(!defined($exon->cdna_start())) {
    throw("exons cdna coding start not defined.");
  }

  if($exon->start_phase() == -1) {
    if($transcript->cdna_coding_start() == $exon->cdna_start()) {
      $exon->start_phase(0);
    }
  } else {
    if($transcript->cdna_coding_start() > $exon->cdna_start()) {
      $exon->start_phase(-1);
    }
  }

  if($exon->end_phase() == -1) {
    if($transcript->cdna_coding_end() == $exon->cdna_end()) {
      # no utr left at end of this exon anymore
      my $cds_len =
        $transcript->cdna_coding_end - $transcript->cdna_coding_start + 1;
      $exon->end_phase($cds_len % 3);
    }
  } else {
    if($transcript->cdna_coding_end() < $exon->cdna_end) {
      # exon end is now utr
      $exon->end_phase(-1);
    }
  }
}



#
# Fixes exons phases of a newly split exon.  The end_phase of the first exon
# and the start_phase of the second exon are set.
#
#
sub set_split_phases {
  my $first_exon = shift;
  my $second_exon = shift;
  my $transcript = shift;

  # need to set first exon end phase and second exon start phase

  if($first_exon->cdna_end > $transcript->cdna_coding_end()) {
    $first_exon->end_phase(-1);
    $second_exon->start_phase(-1);
    return;
  }

  if($first_exon->cdna_end() < $transcript->cdna_coding_start()) {
    $first_exon->end_phase(-1);

    # beginning of CDS could be right at start of second exon
    if($second_exon->cdna_start() == $transcript->cdna_coding_start()) {
      $second_exon->start_phase(0);
    } else {
      $second_exon->start_phase(-1);
    }

    return;
  }

  my $phase;

  if($first_exon->cdna_start() < $transcript->cdna_coding_start()) {
    # first exon is partially 5prime UTR

    my $coding_len = $first_exon->cdna_end() - $transcript->cdna_coding_start();
    $phase = $coding_len % 3;
  } else {
    # first exon should be all CDS

    my $sphase;

    # sometimes start phase may be -1 even though this is all CDS.
    # this is because we have not fixed the start phase yet and it should
    # be 0 due to deletion of the utr
    if($first_exon->start_phase() == -1) {
      $sphase = 0;
    } else {
      $sphase = $first_exon->start_phase();
    }

    $phase = ($first_exon->length() + $sphase) % 3;
  }

  $first_exon->end_phase($phase);
  $second_exon->start_phase($phase);

  return;
}



1;
