
use strict;
use warnings;

package Transcript;

#
# A set of utility methods for dealing with Interim transcripts
# and creating real transcripts out of them.
#


use StatMsg;
use InterimTranscript;
use Utils qw(print_exon);
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Utils::Exception qw(throw info);


use constant MAX_INTRON_LEN => 2e6;


#
# sanity checks the interim exons, and splits this
# interim transcript into parts
#

sub check_iexons {
  my $itranscript = shift;
  my $itranscript_array   = shift;

  my $prev_start = undef;
  my $prev_end = undef;
  my $transcript_seq_region = undef;
  my $transcript_strand     = undef;

  info("checking exons for : " . $itranscript->stable_id());

  my $first = 1;

  foreach my $iexon (@{$itranscript->get_all_Exons}) {

    if ($iexon->fail() || $iexon->is_fatal()) {
      info("  failed/fatal exon, splitting transcript");

      # split this transcript in two, and restart the processing
      # at the beginning of the second transcript
      my $first_transcript = split_itrans($itranscript, $iexon);

      if ($first_transcript) {
        $itranscript_array ||= [];
        push @$itranscript_array, $first_transcript;
      }
      return check_iexons($itranscript, $itranscript_array);
    }

    #sanity check: expect first exon to have cdna_start = 1
    if($first && $iexon->cdna_start != 1) {
      print_exon($iexon);
      throw("Unexpected: first exon does not have cdna_start = 1");
    }
    $first = 0;

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

      info("  scaffold span, splitting transcript");

      my $stat_msg = StatMsg->new(StatMsg::TRANSCRIPT|StatMsg::SCAFFOLD_SPAN);
      $itranscript->add_StatMsg($stat_msg);

      my $keep_exon = 1;
      my $first_transcript = split_itrans($itranscript, $iexon, $keep_exon);

      if ($first_transcript) {
        $itranscript_array ||= [];
        push @$itranscript_array, $first_transcript;
      }
      return check_iexons($itranscript, $itranscript_array);
    }


    # watch out for exons that come in the wrong order

    if((defined($prev_end) && $iexon->strand() == 1 &&
        $prev_end > $iexon->start()) ||
       (defined($prev_start) && $iexon->strand() == -1 &&
        $prev_start < $iexon->end())) {

      info("  inversion, splitting transcript");

      my $stat_msg = StatMsg->new(StatMsg::TRANSCRIPT | StatMsg::INVERT);
      $itranscript->add_StatMsg($stat_msg);

      my $keep_exon = 1;
      my $first_transcript = split_itrans($itranscript, $iexon, $keep_exon);

      if ($first_transcript) {
        $itranscript_array ||= [];
        push @$itranscript_array, $first_transcript;
      }
      return check_iexons($itranscript, $itranscript_array);
    }

    if (!defined($transcript_strand)) {
      $transcript_strand = $iexon->strand();
    } elsif ($transcript_strand != $iexon->strand()) {
      info("  strand flip, splitting transcript");
      my $stat_msg = StatMsg->new(StatMsg::TRANSCRIPT | StatMsg::STRAND_FLIP);
      $itranscript->add_StatMsg($stat_msg);

      my $keep_exon = 1;
      my $first_transcript = split_itrans($itranscript, $iexon, $keep_exon);

      if ($first_transcript) {
        $itranscript_array ||= [];
        push @$itranscript_array, $first_transcript;
      }
      return check_iexons($itranscript, $itranscript_array);
    }

    # watch out for extremely long introns
    my $intron_len = 0;

    if(defined($prev_start)) {
      if($iexon->strand() == 1) {
        $intron_len = $iexon->start - $prev_end + 1;
      } else {
        $intron_len = $prev_start - $iexon->end + 1;
      }
    }

    if($intron_len > MAX_INTRON_LEN) {
      info("  very long intron, splitting transcripts");
      my $keep_exon = 1;
      my $first_transcript = split_itrans($itranscript, $iexon, $keep_exon);

      if($first_transcript) {
        push @$itranscript_array, $first_transcript;
      }
      return check_iexons($itranscript, $itranscript_array);
    }

    $prev_end = $iexon->end();
    $prev_start = $iexon->start();
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
    info("  no exons left in transcript");
  }

  return $itranscript_array;
}


#
# Splits an interim transcript into two parts by discarding
# a 'bad exon' in the middle.
#
# If the 'bad exon' was the first exon then undef is returned, otherwise
# the first transcript is returned.
#
# The passed in transcript is adjusted to become the second transcript
#
# The keep_exon argument is a flag.  If true the 'bad exon' is kept as part
# of the second transcript, otherwise it is discarded.
#
sub split_itrans {
  my $itrans = shift;
  my $bad_exon = shift; #fulcrum exon to split on
  my $keep_exon = shift;

  my @remaining_exons = @{$itrans->get_all_Exons()};
  my @first_exons;
  my $first_trans = InterimTranscript->new();
  $first_trans->stable_id($itrans->stable_id());
  $first_trans->version($itrans->version);

  my $cur_exon = shift(@remaining_exons);

  # all the exons until the 'bad exon' are loaded into the first transcript

  info("==FIRST TRANSCRIPT:\n");
  while($cur_exon && $cur_exon != $bad_exon) {
    print_exon($cur_exon);
    push @first_exons, $cur_exon;
    $first_trans->add_Exon($cur_exon);
    $cur_exon = shift(@remaining_exons);
  }

  if(!$cur_exon) {
    throw("unexpected: could not find bad exon in transcript");
  }

  info("==BAD EXON: ". (($keep_exon) ? 'keeping' : 'discarding'));
  print_exon($bad_exon);

  # keep the 'bad exon' if the flag was set
  if($keep_exon) {
    unshift @remaining_exons, $bad_exon;
  }

  # the remaining exons are reloaded into the second (original) transcript
  $itrans->flush_Exons();
  info("==SECOND TRANSCRIPT:\n");
  foreach my $exon (@remaining_exons) {
    print_exon($exon);
    $itrans->add_Exon($exon);
  }

  #
  # set the coding end/start of the first transcript
  #
  if(@first_exons) {
    my $first_ex = $first_exons[0];
    my $last_ex = $first_exons[$#first_exons];

    # set the coding start and coding end to same as original transcript,
    # then adjust for ommitted exons
    $first_trans->cdna_coding_start($itrans->cdna_coding_start());
    $first_trans->cdna_coding_end($itrans->cdna_coding_end());

    if($last_ex->cdna_end()    < $first_trans->cdna_coding_start() ||
       $first_ex->cdna_start() > $first_trans->cdna_coding_end()) {
      # entire transcript is UTR and is non-coding
      $first_trans->cdna_coding_start(1);
      $first_trans->cdna_coding_end(0);
    }
    elsif($last_ex->cdna_end() >= $first_trans->cdna_coding_start() &&
          $last_ex->cdna_end() < $first_trans->cdna_coding_end()) {
      # coding sequence is cut by coding end
      $first_trans->cdna_coding_end($last_ex->cdna_end());
    }
  }

  #
  # adjust the coding end/start of the second transcript
  #
  if(@remaining_exons) {
    my $first_ex = $remaining_exons[0];
    my $last_ex  = $remaining_exons[$#remaining_exons];

    # all of the cdna coodinates need to be shifted up by the
    # amount of bases used by the first transcript & bad exon
    my $cdna_shift = 0;
    foreach my $ex (@first_exons) {
      $cdna_shift += $ex->length();
      info("cdna_shift = $cdna_shift\n");
    }

    # if we threw away the 'bad' exon we want to take into account
    # how much sequence it had,
    if(!$keep_exon &&
       defined($bad_exon->cdna_start()) &&
       defined($bad_exon->cdna_end())) {
      $cdna_shift += $bad_exon->length();
    }

    info("Shifting CDNA of second transcript by $cdna_shift\n");

    if($cdna_shift) {
      foreach my $ex (@remaining_exons) {
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

  return (@first_exons) ? $first_trans : undef;
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
  $transcript->version($itrans->version);

  info("making final transcript for ". $itrans->stable_id);

  my $translation;

  # the whole translation may have been deleted
  # discard translation if mrna is less than a codon in length
  if($itrans->cdna_coding_end - $itrans->cdna_coding_start + 1 < 3) {
    $translation = undef;
  } else {
    $translation = Bio::EnsEMBL::Translation->new();
    $transcript->translation($translation);
  }

  foreach my $iexon (@{$itrans->get_all_Exons}) {
    my $slice =
      $slice_adaptor->fetch_by_region('chromosome', $iexon->seq_region,
                                     undef, undef,undef, 'CHIMP1');

    my $exon = Bio::EnsEMBL::Exon->new
      (-START     => $iexon->start(),
       -END       => $iexon->end(),
       -STRAND    => $iexon->strand(),
       -PHASE     => $iexon->start_phase(),
       -END_PHASE => $iexon->end_phase(),
       -STABLE_ID => $iexon->stable_id(),
       -SLICE     => $slice);

    $transcript->add_Exon($exon);

    # see if this exon is the start or end exon of the translation
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
        my $translation_end =
          $itrans->cdna_coding_end() - $iexon->cdna_start() + 1;
        $translation->end_Exon($exon);
        $translation->end($translation_end);
      }
    }
  }

  if($translation && !$translation->start_Exon()) {
    print STDERR "Could not find translation start exon in transcript.\n";
    print STDERR "FIRST EXON:\n";
    print_exon($itrans->get_all_Exons->[0]);
    print STDERR "LAST EXON:\n";
    print_exon($itrans->get_all_Exons->[-1], $itrans);
    throw("Unexpected: Could not find translation start exon in transcript\n");
  }
  if($translation && !$translation->end_Exon()) {
    print STDERR "Could not find translation end exon in transcript.\n";
    print STDERR "FIRST EXON:\n";
    print_exon($itrans->get_all_Exons->[0]);
    print STDERR "LAST EXON:\n";
    print_exon($itrans->get_all_Exons->[-1], $itrans);
    throw("Unexpected: Could not find translation end exon in transcript\n");
  }

  return $transcript;
}



1;
