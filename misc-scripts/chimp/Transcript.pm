
use strict;
use warnings;

package Transcript;

#
# A set of utility methods for dealing with Interim transcripts
# and creating real transcripts out of them.
#


use StatMsg;
use InterimTranscript;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Utils::Exception qw(throw);

#
# Sets the phases of the interim exons in an interim transcript
# failed exons are skipped completely.
#
# Note that the phases of the original exon are ignored, so there is a
# possibility of bad translations when the original start exon had a non-zero 
# start phase.
#

sub set_iexon_phases {
  my $itranscript = shift;

  debug("setting exon phases for : " . $itranscript->stable_id());

  my @iexons = @{$itranscript->get_all_Exons()};

  # If all of the CDS was deleted, then all phases are -1

  if ($itranscript->cdna_coding_start == $itranscript->cdna_coding_end + 1) {
    foreach my $ex (@iexons) {
      $ex->start_phase(-1);
      $ex->end_phase(-1);
    }
    return;
  }

  my $cdna_start = 1;
  my $cur_phase  = undef;

  foreach my $iexon (@iexons) {
    next if($iexon->fail());

    my ($start_phase, $end_phase);
    my $cdna_end = $cdna_start + $iexon->length()-1;

    if (defined($cur_phase)) {
      if ($cdna_start <= $itranscript->cdna_coding_end()) {
        $start_phase = $cur_phase; #start phase is last exons end phase
      } else {
        #the end of the last exon was the end of the CDS
        $start_phase = -1;
      }
    } else {
      #sanity check
      if ($cdna_start > $itranscript->cdna_coding_start() &&
          $cdna_start < $itranscript->cdna_coding_end()) {
        throw("Unexpected.  Start of CDS is not in exon?\n" .
              "  exon_cdna_start = $cdna_start\n" .
              "  cdna_coding_start = ".$itranscript->cdna_coding_start());
      }
      if ($cdna_start == $itranscript->cdna_coding_start()) {
        $start_phase = 0;
      } else {
        $start_phase = -1;
      }
    }

    if ($cdna_end < $itranscript->cdna_coding_start() ||
        $cdna_end > $itranscript->cdna_coding_end()) {
      #the end of this exon is outside the CDS
      $end_phase = -1;
    } else {
      #the end of this exon is in the CDS
      #figure out how much coding sequence in the exon
	
      my $coding_start;
      if ($itranscript->cdna_coding_start() > $cdna_start) {
        $coding_start = $itranscript->cdna_coding_start();
      } else {
        $coding_start = $cdna_start;
      }
      my $coding_len = $cdna_end - $coding_start + 1;
	
      if ($start_phase > 0) {
        $coding_len += $start_phase;
      }

      $end_phase = $coding_len % 3;
    }

    $iexon->start_phase($start_phase);
    $iexon->end_phase($end_phase);

    $cdna_start = $cdna_end + 1;
    $cur_phase = ($end_phase >= 0) ? $end_phase : undef;
  }
}


#
# sanity checks the interim exons, and splits this
# interim transcript into parts
#

sub check_iexons {
  my $itranscript = shift;
  my $itranscript_array   = shift;

  my $prev_end = 0;
  my $transcript_seq_region = undef;
  my $transcript_strand     = undef;

  debug("checking exons for : " . $itranscript->stable_id());

  foreach my $iexon (@{$itranscript->get_all_Exons}) {

    #
    # TBD being 'fatal' can have some knockon effects
    # such as the next exon being part of the error (when split due to
    # frameshifting). Do something about this...
    #
    if ($iexon->fail() || $iexon->is_fatal()) {
      debug("  failed/fatal exon, splitting transcript");

      # split this transcript in two, and restart the processing
      # at the beginning of the second transcript
      my $first_transcript = split_itrans($itranscript, $iexon);

      if ($first_transcript) {
        $itranscript_array ||= [];
        push @$itranscript_array, $first_transcript;
      }
      return check_iexons($itranscript, $itranscript_array);
    }

    # sanity check: start must be less than or equal to end
    if ($iexon->end() < $iexon->start()) {
      throw("Unexpected: exon start less than end:\n" .
            $iexon->stable_id().": ".$iexon->start().'-'.$iexon->end());
    }

    # sanity check: cdna length must equal length
    if($iexon->length != $iexon->cdna_end - $iexon->cdna_start + 1) {
      throw("Unexpected: exon cdna length != exon length:\n" .
            $iexon->stable_id().": ".$iexon->start().'-'.$iexon->end() ."\n" .
            "                 " . $iexon->cdna_start.'-'.$iexon->cdna_end());
    }

    if (!defined($transcript_seq_region)) {
      $transcript_seq_region = $iexon->seq_region();
    } elsif ($transcript_seq_region ne $iexon->seq_region()) {

      debug("  scaffold span, splitting transcript");

      my $stat_msg = StatMsg->new(StatMsg::TRANSCRIPT|StatMsg::SCAFFOLD_SPAN);
      $itranscript->add_StatMsg($stat_msg);

      my $first_transcript = split_itrans($itranscript, $iexon);

      if ($first_transcript) {
        $itranscript_array ||= [];
        push @$itranscript_array, $first_transcript;
      }
      return check_iexons($itranscript, $itranscript_array);
    }


    if ($prev_end > $iexon->start()) {
      # something funny has happened, this exon starts before end of previous
      # exon

      debug("  inversion, splitting transcript");

      my $stat_msg = StatMsg->new(StatMsg::TRANSCRIPT | StatMsg::INVERT);
      $itranscript->add_StatMsg($stat_msg);

      my $first_transcript = split_itrans($itranscript, $iexon);

      if ($first_transcript) {
        $itranscript_array ||= [];
        push @$itranscript_array, $first_transcript;
      }
      return check_iexons($itranscript, $itranscript_array);
    }

    if (!defined($transcript_strand)) {
      $transcript_strand = $iexon->strand();
    } elsif ($transcript_strand != $iexon->strand()) {
      debug("  strand flip, splitting transcript");
      my $stat_msg = StatMsg->new(StatMsg::TRANSCRIPT | StatMsg::STRAND_FLIP);
      $itranscript->add_StatMsg($stat_msg);

      my $first_transcript = split_itrans($itranscript, $iexon);

      if ($first_transcript) {
        $itranscript_array ||= [];
        push @$itranscript_array, $first_transcript;
      }
      return check_iexons($itranscript, $itranscript_array);
    }
  }

  $itranscript_array ||= [];

  #
  # if there exons left after all the splitting,
  # then add this transcript to the array
  #
  my $total_exons = scalar(@{$itranscript->get_all_Exons});
  if ($total_exons > 0) {
    push @$itranscript_array, $itranscript;
  } else {
    debug("  no exons left in transcript");
  }

  return $itranscript_array;
}


#
# splits an interim transcript into two parts by discarding
# a 'bad exon' in the middle.
# If the 'bad exon' was the first exon then undef is returned, otherwise
# the first transcript is returned.
# the passed in transcript is adjusted to become the second transcript
#
sub split_itrans {
  my $itrans = shift;
  my $bad_exon = shift; #fulcrum exon to split on

  my @all_exons = @{$itrans->get_all_Exons()};
  my @first_exons;
  my @second_exons;
  my $first_trans = InterimTranscript->new();
  $first_trans->stable_id($itrans->stable_id());

  my $cur_exon = shift(@all_exons);

  # all the exons until the 'bad exon' are loaded into the first transcript

  while($cur_exon && $cur_exon != $bad_exon) {
    push @first_exons, $cur_exon;
    $first_trans->add_Exon($cur_exon);
    $cur_exon = shift(@all_exons);
  }

  if(!$cur_exon) {
    throw("unexpected: could not find bad exon in transcript");
  }

  # the remaining exons are reloaded into the second (original) transcript
  $itrans->flush_Exons();
  while(my $ex = shift(@all_exons)) {
    push @second_exons, $ex;
    $itrans->add_Exon($ex);
  }

  #
  # set the coding end/start of the first transcript
  #
  if(@first_exons) {
    my $first_ex = $first_exons[0];
    my $last_ex = $first_exons[$#first_exons];

    # set the coding start and coding end to same as original transcript,
    # then adjust for ommitted
    $first_trans->cdna_coding_start($itrans->cdna_coding_start());
    $first_trans->cdna_coding_end($itrans->cdna_coding_end());

    if($last_ex->cdna_end()    < $first_trans->cdna_coding_start() ||
       $first_ex->cdna_start() > $first_trans->cdna_coding_end()) {
      # entire transcript is UTR and is non-coding
      $first_trans->cdna_coding_start(1);
      $first_trans->cdna_coding_end(0);
    }
    elsif($last_ex->cdna_end() > $first_trans->cdna_coding_start() &&
	  $last_ex->cdna_end() < $first_trans->cdna_coding_end()) {
      # coding sequence is cut by coding end
      $first_trans->cdna_coding_end($last_ex->cdna_end());
    }
  }

  #
  # adjust the coding end/start of the second transcript
  #
  if(@second_exons) {
    my $first_ex = $second_exons[0];
    my $last_ex  = $second_exons[$#second_exons];

    # all of the cdna coodinates need to be shifted up by the
    # amount of bases used by the first transcript & bad exon
    my $cdna_shift = 0;
    foreach my $ex (@first_exons) {
      $cdna_shift += $ex->length();
    }
    if(!$bad_exon->fail()) {
      $cdna_shift += $bad_exon->length();
    }

    if($cdna_shift) {
      foreach my $ex (@second_exons) {
        if(!$ex->fail()) {
          $ex->cdna_start($ex->cdna_start() - $cdna_shift);
          $ex->cdna_end($ex->cdna_end() - $cdna_shift);
        }
      }
      $itrans->move_cdna_coding_start(-$cdna_shift);
      $itrans->move_cdna_coding_end(-$cdna_shift);

      if($itrans->cdna_coding_start() < 1) {
        # transcript was split in middle of cds
        $itrans->cdna_coding_start(1);
      }
      if($itrans->cdna_coding_end() < 1) {
        #there is no coding sequence left in this transcript
        $itrans->cdna_coding_end(0);
      }
    }
  }

  return $first_trans if(@first_exons);

  return undef;
}

#
# creates proper ensembl transcripts and exons from interim transcripts
# and exons.
#
sub make_Transcript {
  my $itrans = shift;
  my $slice_adaptor = shift;

  my $transcript = Bio::EnsEMBL::Transcript->new();
  $transcript->stable_id($itrans->stable_id);

  debug("making final transcript for ". $itrans->stable_id);

  my $translation;

  # the whole translation may have been deleted
  if ($itrans->cdna_coding_start == $itrans->cdna_coding_end + 1) {
    $translation = undef;
  } else {
    $translation = Bio::EnsEMBL::Translation->new();
    $transcript->translation($translation);
  }

  foreach my $iexon (@{$itrans->get_all_Exons}) {
    my $slice =
      $slice_adaptor->fetch_by_region('scaffold', $iexon->seq_region);

    my $exon = Bio::EnsEMBL::Exon->new
      (-START     => $iexon->start(),
       -END       => $iexon->end(),
       -STRAND    => $iexon->strand(),
       -PHASE     => $iexon->start_phase(),
       -END_PHASE => $iexon->end_phase(),
       -STABLE_ID => $iexon->stable_id(),
       -SLICE     => $slice);

    $transcript->add_Exon($exon);

    #
    # see if this exon is the start or end exon of the translation
    #

    if ($translation) {
      if ($iexon->cdna_start() <= $itrans->cdna_coding_start() &&
          $iexon->cdna_end()   >= $itrans->cdna_coding_start()) {
        my $translation_start =
          $itrans->cdna_coding_start() - $iexon->cdna_start() + 1;
        $translation->start_Exon($exon);
        $translation->start($translation_start);
      }

      if ($iexon->cdna_start() <= $itrans->cdna_coding_end() &&
          $iexon->cdna_end()   >= $itrans->cdna_coding_end()) {
        debug("end exon=".$exon->stable_id()."\n".
              "  ex_coding_start=".$iexon->cdna_start()."\n".
              "  ex_coding_end=".$iexon->cdna_end()."\n".
              "  start=".$iexon->start()."\n".
              "  end=".$iexon->end()."\n".
              "  transl_end = ".$itrans->cdna_coding_end());

        my $translation_end =
          $itrans->cdna_coding_end() - $iexon->cdna_start() + 1;
        $translation->end_Exon($exon);
        $translation->end($translation_end);
      }
    }
  }

  return $transcript;
}


sub debug {
  my $msg = shift;
  print STDERR $msg, "\n";
}


1;
