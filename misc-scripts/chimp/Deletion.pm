use strict;
use warnings;

package Deletion;

use InterimExon;
use StatMsg;
use Length;

use Bio::EnsEMBL::Utils::Exception qw(throw info);


###############################################################################
# process_delete
#
# processes a deletion in an exon
###############################################################################

sub process_delete {
  my $cdna_del_pos_ref  = shift;
  my $del_len           = shift;
  my $exon        = shift;
  my $transcript  = shift;
  my $entire_delete = shift;

  my $del_start = $$cdna_del_pos_ref + 1;
  my $del_end   = $del_start + $del_len - 1;

  info((($entire_delete) ? 'entire ' : '')."delete ($del_len)");
  #print STDERR "BEFORE cds: ", $transcript->cdna_coding_start, '-', 
  #             $transcript->cdna_coding_end(), "\n";
  #print STDERR "BEFORE del_start = $del_start\n";

  # sanity check, deletion should be completely in
  # or adjacent to exon boundaries
  if(!$entire_delete && ($del_start < $exon->cdna_start() - 1 ||
			 $del_start > $exon->cdna_end() + 1)) {

    throw("Unexpected: deletion is outside of exon boundary\n" .
          "     del_start       = $del_start\n" .
          "     cdna_exon_start =". $exon->cdna_start() .
          "     cdna_exon_end   =". $exon->cdna_end());
  }

  # break delete into composite parts and deal with each part seperately

  #
  # deal with five prime UTR portion of delete
  #
  if($del_start < $transcript->cdna_coding_start()) {
    my $utr_del_len;

    if($del_end >= $transcript->cdna_coding_start()) {
      $utr_del_len = $transcript->cdna_coding_start() - $del_start;
    } else {
      $utr_del_len = $del_len;
    }

    process_five_prime_utr_delete($cdna_del_pos_ref, $utr_del_len,
				  $exon, $transcript);

    # take away the processed part of the deletion
    $del_start = $$cdna_del_pos_ref + 1;
    $del_len   -= $utr_del_len;
    $del_end   = $del_start + $del_len - 1;
  }

  return if($del_len == 0); # no deletion left

  #
  # deal with CDS portion of delete
  #
  if($del_end >= $transcript->cdna_coding_start() &&
     $del_start <= $transcript->cdna_coding_end()) {

    my $cds_del_len;

    if($del_end > $transcript->cdna_coding_end()) {
      $cds_del_len = $transcript->cdna_coding_end() - $del_start + 1;
    } else {
      $cds_del_len = $del_len;
    }

    process_cds_delete($cdna_del_pos_ref, $cds_del_len, $exon, $transcript,
		      $entire_delete);

    # take away the processed part of the deletion
    # the cdna start is in the same place because
    $del_start += $$cdna_del_pos_ref + 1;
    $del_len  -= $cds_del_len;
    $del_end   = $del_start + $del_len - 1;
  }

  return if($del_len == 0); # no deletion left

  #
  # deal with 3prime portion of delete
  #

  # sanity check:
  if($del_start <= $transcript->cdna_coding_end()) {
    throw("Unexpected. 3' UTR delete starts before coding end.");
  }

  process_three_prime_utr_delete($cdna_del_pos_ref, $del_len, $exon,
				 $transcript);

  return;
}

###############################################################################
# process_five_prime_utr_delete
#
# processes a deletion in the five prime utr of a transcript
###############################################################################

sub process_five_prime_utr_delete {
  my $cdna_del_pos_ref = shift;
  my $del_len          = shift;
  my $exon             = shift;
  my $transcript       = shift;

  info("delete ($del_len) in 5' utr");

  # shift up the CDS
  $transcript->move_cdna_coding_start(-$del_len);
  $transcript->move_cdna_coding_end(-$del_len);

  # create a status message and add it to the exon
  my $code = StatMsg::EXON | StatMsg::DELETE | StatMsg::FIVE_PRIME |
             StatMsg::UTR | Length::length2code($del_len);
  $exon->add_StatMsg(StatMsg->new($code));

  return;
}

###############################################################################
# process_three_prime_utr_delete
#
# processes a deletion in the three prime utr of a transcript
###############################################################################

sub process_three_prime_utr_delete {
  my $cdna_del_pos_ref = shift;
  my $del_len          = shift;
  my $exon             = shift;
  my $transcript       = shift;

  #do not have to do anything...
  info("delete ($del_len) in 3' utr");

  # create a status message and add it to the exon
  my $code = StatMsg::EXON | StatMsg::DELETE | StatMsg::THREE_PRIME |
             StatMsg::UTR | Length::length2code($del_len);
  $exon->add_StatMsg(StatMsg->new($code));

  return;
}

###############################################################################
# process_cds_delete
#
# processes a deletion in the cds of a transcript
###############################################################################

sub process_cds_delete {
  my $cdna_del_pos_ref = shift;
  my $del_len          = shift;
  my $exon             = shift;
  my $transcript       = shift;
  my $entire_delete    = shift;

  info("delete ($del_len) in cds");

  my $del_start = $$cdna_del_pos_ref + 1;
  my $del_end   = $del_start + $del_len - 1;

  my $code = StatMsg::EXON | StatMsg::DELETE | StatMsg::CDS |
             Length::length2code($del_len);

  my $frameshift = $del_len % 3;


  #
  # case 1: delete is all of CDS
  #
  if($del_start == $transcript->cdna_coding_start() &&
     $del_end   == $transcript->cdna_coding_end()) {
    info("delete ($del_len) is all of cds");

    $code |= StatMsg::ENTIRE;

    # move up CDS end to account for CDS deletion
    $transcript->move_cdna_coding_end(-$del_len);
  }

  #
  # case 2: delete is at start of CDS
  #
  elsif($del_start == $transcript->cdna_coding_start()) {
    info("delete ($del_len) at start of cds");

    $code |= StatMsg::FIVE_PRIME;

    # move up CDS end to account for CDS deletion
    $transcript->move_cdna_coding_end(-$del_len);

    if($frameshift) {
      $code |= StatMsg::FRAMESHIFT if($frameshift);

      # move down CDS start to put reading frame back (shrink CDS)
      info("shifting cds start to restore reading frame");
      $transcript->move_cdna_coding_start(3 - $frameshift);
    }
  }

  #
  # case 3: delete is at end of CDS
  #
  elsif($del_end == $transcript->cdna_coding_end()) {
    info("delete ($del_len) at end of cds");

    $code |= StatMsg::THREE_PRIME;

    # move up CDS end to account for CDS deletion
    $transcript->move_cdna_coding_end(-$del_len);

    if($frameshift) {
      $code |= StatMsg::FRAMESHIFT if($frameshift);

      # move up CDS end to put reading frame back (shrink CDS)
      info("shifting cds end to restore reading frame");
      $transcript->move_cdna_coding_end(3 - $frameshift);
    }
  }

  #
  # case 4: delete is in middle of CDS
  #
  elsif($del_end   > $transcript->cdna_coding_start() &&
        $del_start < $transcript->cdna_coding_end()) {
    info("delete ($del_len) in middle of cds");

    $code |= StatMsg::MIDDLE;

    # move up CDS end to account for CDS deletion
    $transcript->move_cdna_coding_end(-$del_len);

    if($frameshift && !$entire_delete) {
      $code |= StatMsg::FRAMESHIFT if($frameshift);

      # this is going to require splitting the exon
      # to make a frameshift deletion

      #first exon is going to end right before deletion
      my $first_len  = $del_start - $exon->cdna_start();
      my $intron_len = 3 - $frameshift;

      #reduce the length of the CDS by the length of the new intron
      $transcript->move_cdna_coding_end(-$intron_len);

      # the next match that is added to the cdna position will have too much
      # sequence because we used part of the sequence to create the frameshift
      # intron, compensate by reducing cdna position by intron len
      $$cdna_del_pos_ref -= $intron_len;

      info("introducing frameshift intron ($intron_len) " .
            "to maintain reading frame");

      # very short exons can be entirely consumed by the intron 
      if($intron_len == $exon->length()) {
        $code |= StatMsg::ALL_INTRON;
        $exon->fail(1);
      } 
      elsif($intron_len > $exon->length()) {
        $code |= StatMsg::CONFUSED | StatMsg::ALL_INTRON;
        $exon->fail(1);
      }
      elsif($first_len + $intron_len >= $exon->length()) {
        # we may have encountered a delete at the very end of the exon
        # in this case we have to take the intron out of the end of this exon
        # since we are not creating a second one

        if($exon->strand() == 1) {
          $exon->end($exon->end - $intron_len);
        } else {
          $exon->start($exon->start + $intron_len);
        }
        $exon->cdna_end($exon->cdna_end - $intron_len);
      } else {
        # second exon is going to start right after 'frameshift intron'

        if($exon->strand == 1) {
          # end the current exon at the beginning of the deletion
          # watch out though, because we may be at the very beginning of
          # the exon in which case we do not want to create one

          if($first_len) {
            my $first_exon = InterimExon->new();

            # copy the original exon and adjust the coords as necessary
            %{$first_exon} = %{$exon};
            $first_exon->cdna_end($exon->cdna_start() + $first_len - 1);
            $first_exon->end($first_exon->start() + $first_len - 1);
            $transcript->add_Exon($first_exon);

            $exon->cdna_start($first_exon->cdna_end() + 1);
            $exon->flush_StatMsgs();
          }

          # start next exon after new intron
          $exon->start($exon->start() + $first_len + $intron_len);
          $exon->cdna_end($exon->cdna_end - $intron_len);
        } else {
          if($first_len) {
            my $first_exon = InterimExon->new();

            # copy the original exon and adjust the coords as necessary
            %{$first_exon} = %{$exon};
            $first_exon->cdna_end($exon->cdna_start + $first_len - 1);
            $first_exon->start($exon->end() - $first_len + 1);
            $transcript->add_Exon($first_exon);

            $exon->cdna_start($first_exon->cdna_end() + 1);
            $exon->flush_StatMsgs();
          }

          # start next exon after new intron
          $exon->end($exon->end() - ($first_len + $intron_len));
          $exon->cdna_end($exon->cdna_end - $intron_len);
        }
      }
    }
  }

  # sanity check:
  else {
    throw("Unexpected: CDS delete appears to be outside of CDS:\n" .
         "  del_start = $del_start\n".
         "  del_end   = $del_end\n" .
         "  cdna_coding_start = ".$transcript->cdna_coding_start() . "\n".
         "  cdna_coding_end   = ".$transcript->cdna_coding_end() . "\n");
  }

  $exon->add_StatMsg(StatMsg->new($code));

  return;
}



#sub print_exon {
  #my $exon = shift;
#
  #print "EXON:\n";
  #print "cdna_start = ". $exon->cdna_start() . "\n";
  #print "cdna_end   = ". $exon->cdna_end() . "\n";
  #print "start             = ". $exon->start() . "\n";
  #print "end               = ". $exon->end() . "\n";
  #print "strand            = ". $exon->strand() . "\n\n";

  #return;
#}

1;
