use strict;
use warnings;

use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;

use InterimTranscript;
use InterimExon;
use StatMsg;
use Deletion;

use Bio::EnsEMBL::Utils::Exception qw(throw);



my $verbose;

{  #block to avoid namespace pollution

  my ($hhost, $hdbname, $huser, $hpass, $hport, $hassembly,  #human vars
      $hchromosome, $hstart, $hend,
      $chost, $cdbname, $cuser, $cpass, $cport, $cassembly,  #chimp vars
      $help);

  GetOptions('hhost=s'   => \$hhost,
             'hdbname=s' => \$hdbname,
             'huser=s'   => \$huser,
             'hpass=s'   => \$hpass,
             'hport=i'   => \$hport,
             'hassembly=s' => \$hassembly,
             'hchromosome=s' => \$hchromosome,
             'hstart=i'  => \$hstart,
             'hend=i'    => \$hend,
             'chost=s'   => \$chost,
             'cdbname=s' => \$cdbname,
             'cuser=s'   => \$cuser,
             'cpass=s'   => \$cpass,
             'cport=i'   => \$cport,
             'cassembly=s' => \$cassembly,
             'help'      => \$help,
             'verbose'   => \$verbose);


  usage() if($help);
  usage("-hdbname option is required") if (!$hdbname);
  usage("-cdbname option is required") if (!$cdbname);

  $hport ||= 3306;
  $cport ||= 3306;

  $hdbname ||= 'localhost';
  $cdbname ||= 'localhost';

  $hassembly ||= 'NCBI34';
  $cassembly ||= 'BROAD1';

  debug("Connecting to human database");

  my $human_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (-host    => $hhost,
     -dbname => $hdbname,
     -pass   => $hpass,
     -user   => $huser,
     -port   => $hport);

  debug("Connecting to chimp database");

  my $chimp_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (-host    => $chost,
     -dbname  => $cdbname,
     -pass    => $cpass,
     -user    => $cuser,
     -port    => $cport);


  $human_db->dnadb($chimp_db);

  my $slice_adaptor   = $human_db->get_SliceAdaptor();
  my $gene_adaptor    = $human_db->get_GeneAdaptor();


  debug("Fetching chromosomes");

  my $slices;

  if($hchromosome) {
    my $slice = $slice_adaptor->fetch_by_region('chromosome',
                                                $hchromosome,
                                                $hstart, $hend, undef,
                                                $hassembly);
    if(!$slice) {
      throw("unknown chromosome $hchromosome");
    }

    $slices = [$slice];
  } else {
    $slices = $slice_adaptor->fetch_all('chromosome', $hassembly);
  }

  my %err_count;
  my %err_descs;
  my %ok_err_count;

  my $total_transcripts = 0;
  my $mapped_transcripts = 0;
  my $total_known   = 0;
  my $mapped_known  = 0;

  foreach my $slice (@$slices) {

    debug("Chromosome: " . $slice->seq_region_name());
    debug("Fetching Genes");

    my $genes = $gene_adaptor->fetch_all_by_Slice($slice);

    foreach my $gene (reverse @$genes) {
      debug("Gene: ".$gene->stable_id);
      my $transcripts = $gene->get_all_Transcripts();

      foreach my $transcript (@$transcripts) {
        next if(!$transcript->translation); #skip pseudo genes

        $total_transcripts++;
        $total_known++ if($transcript->is_known);

        my @mapped_trans = transfer_transcript($human_db,$chimp_db,$transcript,
                                             $hassembly, $cassembly);

        if(@mapped_trans) {
          foreach my $mapped_trans (@mapped_trans) {
            debug("\nTranslation:" . $mapped_trans->translate()->seq()."\n\n");
          }
          $mapped_known++ if($transcript->is_known());
          $mapped_transcripts++;

          while(my ($ec, $edesc) = pop_err()) {
            $ok_err_count{$ec}++;
          }

        } else {
          debug("Transcript failed to map\n");
          while(my ($ec, $edesc) = pop_err()) {
            next if($ec == ErrCode::SHORT_FRAMESHIFT_DELETE);
            next if($ec == ErrCode::SHORT_FRAMESHIFT_INSERT);
            $err_count{$ec}++;
            if($edesc) {
              $err_descs{$ec} ||= [];
              push @{$err_descs{$ec}}, $edesc;
            }
          }
        }

        print STDERR "Total Transcripts : $total_transcripts\n";
        print STDERR "Known  Transcripts: $total_known\n";
        print STDERR "Mapped Transcripts: $mapped_transcripts\n";
        print STDERR "Known Mapped Transcripts: $mapped_known\n";

      }
    }
  }

  print STDERR "Mapping Failure Summary\n";
  print STDERR "#Occurances\tDescription\n";
  print STDERR "-----------\t----------------\n";

  foreach my $ec (sort {$err_count{$a} <=> $err_count{$b}} keys %err_count) {
    print STDERR $err_count{$ec}, "\t\t", ec2str($ec), "\n";
    #foreach my $desc (@{$err_descs{$ec}}) {
    #  print STDERR "\t\t  ", $desc, "\n";
    #}
  }

  print STDERR "Other Info\n";
    foreach my $ec (sort {$ok_err_count{$a} <=> $ok_err_count{$b}} 
                    keys %ok_err_count) {
    print STDERR $ok_err_count{$ec}, "\t\t", ec2str($ec), "\n";
  }

}


###############################################################################
#
# transfer_transcript
#
###############################################################################

sub transfer_transcript {
  my $human_db   = shift;
  my $chimp_db   = shift;
  my $transcript = shift;
  my $hassembly  = shift;
  my $cassembly  = shift;

  debug("Transcript: " . $transcript->stable_id());

  my $cs_adaptor      = $human_db->get_CoordSystemAdaptor();
  my $asmap_adaptor   = $human_db->get_AssemblyMapperAdaptor();

  my $chimp_cs = $cs_adaptor->fetch_by_name('scaffold',  $cassembly);
  my $chimp_ctg_cs = $cs_adaptor->fetch_by_name('contig');
  my $human_cs = $cs_adaptor->fetch_by_name('chromosome', $hassembly);

  my $mapper = $asmap_adaptor->fetch_by_CoordSystems($chimp_cs, $human_cs);
  my $scaf2ctg_mapper = $asmap_adaptor->fetch_by_CoordSystems($chimp_cs,
                                                              $chimp_ctg_cs);

  my $human_exons = $transcript->get_all_Exons();

  if(!$transcript->translation()) { #watch out for pseudogenes
    debug("pseudogene - discarding");
    return ();
  }

  my $trans_mapper = $transcript->get_TranscriptMapper();

  my $chimp_cdna_pos = 0;
  my $cdna_exon_start = 1;

  my $chimp_transcript = InterimTranscript->new();
  $chimp_transcript->stable_id($transcript->stable_id());
  $chimp_transcript->cdna_coding_start($transcript->cdna_coding_start());
  $chimp_transcript->cdna_coding_end($transcript->cdna_coding_end());

  my @chimp_exons;

 EXON:
  foreach my $human_exon (@$human_exons) {

    my $chimp_exon = InterimExon->new();
    $chimp_exon->stable_id($stable_id);
    $chimp_exon->cdna_start($cdna_exon_start);

    my @coords = $mapper->map($exon->seq_region_name, $exon->seq_region_start,
                              $exon->seq_region_end, $exon->seq_region_strand,
                              $human_cs);

    if(@coords == 1) {
      my $c = $coords[0];

      if($c->isa('Bio::EnsEMBL::Mapper::Gap')) {
        #
        # complete failure to map exon
        #
	my $stat_code = StatMsg::EXON | StatMsg::DELETE | StatMsg::ENTIRE;

        $chimp_exon->fail(1);
	$chimp_transcript->add_Exon($chimp_exon);

        #if this is a UTR exon, this is ok
        my @cds_coords = $trans_mapper->genomic2cds($c->start(),$c->end(),
                                                    $exon->strand);

        #check if this exon was entirely UTR
        if(@cds_coords == 1 &&
           $cds_coords[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
          debug('entire utr exon deletion - ok');

	  $stat_code |= StatMsg::UTR;

	  $chimp_exon->add_StatMsg(StatMsg->new($stat_code));

          # we can simply throw out this exon, it was all UTR
        } else {
          ### TBD: can do more sensible checking and maybe split transcript
	  $stat_code |= StatMsg::CDS;

	  if(@cds_coords > 1) {
	    # some UTR was deleted too
	    $stat_code |= StatMsg::UTR;
	  }

          $chimp_exon->add_StatMsg(StatMsg->new($stat_code));

          return ();
        }
      } else {
        #
        # exon mapped completely
        #

        $chimp_exon->start($c->start());
        $chimp_exon->end($c->end());
        $chimp_exon->cdna_start($cdna_exon_start);
        $chimp_exon->cdna_end($cdna_exon_start + $chimp_exon->length() - 1);
        $chimp_exon->strand($c->strand());
        $chimp_exon->seq_region($c->id());

        $chimp_cdna_pos += $c->length();
        $cdna_exon_start += $c->length();

        $chimp_transcript->add_Exon($chimp_exon);
      }
    } else {
      my $num = scalar(@coords);

      get_coords_extent(\@coords, $chimp_exon);

      if($chimp_exon->fail()) {
        #failed to obtain extent of coords due to scaffold spanning
        #strand flipping, or exon inversion

	$transcript->add_Exon($chimp_exon);

        ### TBD - fix this, may be ok to drop exon and continue esp. if exon
        ###       is entirely UTR
        return ();
      }

      for(my $i=0; $i < $num; $i++) {
        my $c = $coords[$i];

        if($c->isa('Bio::EnsEMBL::Mapper::Gap')) {

          #
          # deletion in chimp, insert in human
          #
          Deletion::process_deletion(\$chimp_cdna_pos, $c->length(),
				     $chimp_exon, $chimp_transcript);

        } else {
          # can actually end up with adjacent inserts and deletions so we need
          # to take the previous coordinate skipping over gaps that may have
          # been first

          my $prev_c = undef;

          for(my $j = $i-1; $j >= 0 && !defined($prev_c); $j--) {
            if($coords[$j]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
              $prev_c = $coords[$j];
            }
          }

          if($prev_c) {

            my $insert_len;
            if($exon_strand == 1) {
              $insert_len = $c->start() - $prev_c->end() - 1;
            } else {
              $insert_len = $prev_c->start() - $c->end() - 1;
            }

            #sanity check:
            if($insert_len < 0) {
              throw("Unexpected - negative insert " .
                    "- undetected exon inversion?");
            }

            if($insert_len > 0) {

              #
              # insert in chimp, deletion in human
              #

              process_insertion(\$chimp_cdna_pos, $insert_len,
                                $chimp_exon, $chimp_transcript);

              $chimp_cdna_pos += $insert_len;
            }
          }

          $chimp_cdna_pos += $c->length();
        }
      }  # foreach coord

      $cdna_exon_start += $chimp_exon->length();

      $chimp_transcript->add_Exon($chimp_exon);
    }
  } # foreach exon

  my $slice_adaptor = $chimp_db->get_SliceAdaptor();

  return create_transcripts(\@chimp_exons,
                            $cdna_coding_start, $cdna_coding_end,
                            $slice_adaptor);
}


###############################################################################
# process_insertion
#
###############################################################################

sub process_insertion {
  my $cdna_ins_pos_ref = shift;   #basepair to left of insert
  my $insert_len       = shift;
  my $exon             = shift;
  my $transcript       = shift;

  # sanity check, insert should be completely in exon boundaries
  if($$cdna_ins_pos_ref <  $exon->cdna_start() ||
     $$cdna_ins_pos_ref >= $exon->cdna_end()) {

    # because some small (<3bp) matches can be completely eaten away by the
    # introduction of frameshift introns it is possible to get an insert
    # immediately before a newly created (i.e.) split intron

    if($$cdna_ins_pos_ref < $exon->cdna_start  &&
       $$cdna_ins_pos_ref + 3 >= $exon->cdna_start  ) {
      ### TBD not sure what should be done with this situation

      $exon->add_StatMsg(StatMsg->new(StatMsg::PART_EXON_CONFUSED));
      $exon->fail(1);
      $transcript->fail(1);
      return;
    }

    throw("Unexpected: insertion is outside of exon boundary\n" .
          "     ins_left       = $$cdna_ins_pos_ref\n" .
          "     ins_right      = " . ($$cdna_ins_pos_ref+1) . "\n" .
          "     cdna_exon_start = ". $exon->cdna_start()."\n" .
          "     cdna_exon_end   = ". $exon->cdna_end()."\n";);
  }


  #
  # case 1: insert in CDS
  #
  if($$cdna_ins_pos_ref >= $transcript->cdna_coding_start() &&
     $$cdna_ins_pos_ref <  $transcript->cdna_coding_end()) {

    debug("insertion in cds ($insert_len)");

    if($insert_len > $MAX_CODING_INDEL) {
      # too much coding sequence has been inserted

      my $sm = StatMsg->new(StatMsg::PART_EXON_CDS_INSERT_TOO_LONG);
      $exon->add_StatMsg($sm);
      $exon->fail(1);
      $transcript->fail(1);
      return;
    }

    # adjust CDS end accordingly
    $transcript->move_cdna_coding_end($insert_len);

    my $frameshift = $insert_len % 3;

    if($frameshift) {
      if($insert_len > $MAX_FRAMESHIFT_INDEL) {
        my $sm=StatMsg->new(StatMsg::PART_EXON_CDS_FRAMESHIFT_INSERT_TOO_LONG);
        $exon->add_StatMsg($sm);
        $exon->fail(1);
        $transcript->fail(1);
        return;
      }

      $exon->add_StatMsg(StatMsg->new(StatMsg::SHORT_FRAMESHIFT_INSERT));

      # need to create frameshift intron to get reading frame back on track
      # exon needs to be split into two

      debug("introducing frameshift intron to maintain reading frame");

      # first exon ends right before insert
      my $first_len  = $$cdna_ins_pos_ref - $exon->cdna_start() + 1;

      # copy the original exon and adjust coords of each to perform 'split'
      my $first_exon = InterimExon->new();
      %{$first_exon} = %{$exon};
      $exon->flush_StatMsgs();

      # frame shift intron eats into start of inserted region
      # second exon is going to start right after 'frameshift intron'
      # which in cdna coords is immediately after last exon
      $first_exon->cdna_end($first_exon->cdna_start + $first_len - 1);
      $exon->cdna_start($first_exon->cdna_end + 1);

      # decrease the length of the CDS by the length of new intron
      $transcript->move_cdna_coding_end(-$frameshift);

      # the insert length will be added to the cdna_position
      # but part of the insert was used to create the intron and is not cdna
      # anymore, so adjust the cdna_position to compensate
      $$cdna_ins_pos_ref -= $frameshift;

      ### TBD may have to check we have not run up to end of CDS here

      if($exon->strand() == 1) {
        # end the first exon at the beginning of the insert
        $first_exon->end($first_exon->start() + $first_len -1 );

        # start the next exon after the frameshift intron
        $exon->start($exon->start() + $first_len + $frameshift);
      } else {
        $first_exon->start($first_exon->end() - $first_len + 1);

        # start the next exon after the frameshift intron
        $exon->end($exon->end() - ($first_len + $frameshift));
      }
      return;
    }
  }

  #
  # case 2: insert in 5 prime UTR (or between 5prime UTR and CDS)
  #
  elsif($$cdna_ins_pos_ref < $$cdna_coding_start_ref) {
    debug("insertion ($insert_len) in 5' utr");

    #shift the coding region down as result of insert
    $transcript->move_cdna_coding_start($insert_len);
    $transcript->move_cdna_coding_end($insert_len);
  }

  #
  # case 3: insert in 3 prime UTR (or between 3prime UTR and CDS)
  #
  elsif($$cdna_ins_pos_ref >= $$cdna_coding_end_ref) {
    debug("insert ($insert_len) in 3' utr");

    #do not have to do anything
  }

  #
  # default: sanity check
  #
  else {
    throw("Unexpected insert case encountered");
  }

  return [];
}

###############################################################################
# get_coords_extent
#
# given a list of coords returns the start, end, strand, seq_region
# of the span of the coords
#
# undef is returned if the coords flip strands, have an inversion,
# or cross multiple seq_regions
#
###############################################################################

sub get_coords_extent {
  my $coords = shift;
  my $chimp_exon = shift;

  my($start, $end, $strand, $seq_region);

  my $stat_code = StatMsg::EXON;

  foreach my $c (@$coords) {
    next if($c->isa('Bio::EnsEMBL::Mapper::Gap'));

    if(!defined($seq_region)) {
      $seq_region = $c->id();
    }
    elsif($seq_region ne $c->id()) {
      $chimp_exon->fail(1);
      $stat_code |= StatMsg::SCAFFOLD_SPAN;
      $chimp_exon->add_StatMsg(StatMsg->new($stat_code));
      return;
    }

    if(!defined($strand)) {
      $strand = $c->strand();
    }
    elsif($strand != $c->strand()) {
      $chimp_exon->fail(1);
      $stat_code |= StatMsg::STRAND_FLIP;
      $chimp_exon->add_StatMsg(StatMsg->new($stat_code));
      return;
    }

    if(!defined($start)) {
      $start = $c->start if(!defined($start));
    } else {
      if($strand == 1 && $start > $c->start()) {
        $chimp_exon->fail(1);
	$stat_code |= StatMsg::INVERT;
        $chimp_exon->add_StatMsg(StatMsg->new($stat_code));
        return;
      }
      if($strand == -1 && $start < $c->start()) {
        $chimp_exon->fail(1);
	$stat_code |= StatMsg::INVERT;
        $chimp_exon->add_StatMsg(StatMsg->new($stat_code));
        return;
      }

      if($start > $c->start()) {
        $start = $c->start();
      }
    }
	
    if(!defined($end)) {
      $end = $c->end();
    } else {
      if($strand == 1 && $end > $c->end()) {
        $chimp_exon->fail(1);
	$stat_code |= StatMsg::INVERT;
        $chimp_exon->add_StatMsg(StatMsg->new($stat_code));
        return;
      }
      if($strand == -1 && $end < $c->end()) {
        $chimp_exon->fail(1);
	$stat_code |= StatMsg::INVERT;
        $chimp_exon->add_StatMsg(StatMsg->new($stat_code));
        return;
      }
      if($c->end > $end) {
        $end = $c->end();
      }
    }
  }

  $chimp_exon->start($start);
  $chimp_exon->end($end);
  $chimp_exon->strand($strand);
  $chimp_exon->seq_region($seq_region);
  $chimp_exon->cdna_end($chimp_exon->cdna_start() + $chimp_exon->length() - 1);
}

################################################################################
# create_transcripts
#
################################################################################

sub create_transcripts {
  my $exon_coords       = shift;
  my $cdna_coding_start = shift;
  my $cdna_coding_end   = shift;
  my $slice_adaptor     = shift;

  my $transcript_strand = undef;
  my $transcript_seq_region = undef;

  my $transcript = Bio::EnsEMBL::Transcript->new();
  my $translation = Bio::EnsEMBL::Translation->new();

  my @chimp_transcripts;

  my $cdna_start = 1;
  my $cur_phase  = undef;

  foreach my $exon_coord (@$exon_coords) {
    my($exon_start, $exon_end, $exon_strand, $seq_region) = @$exon_coord;

    #sanity check:
    if($exon_end < $exon_start) {
      my $exon_dump = '';
      foreach my $ex (@$exon_coords) {
        if($ex->[0] > $ex->[1]) {
          $exon_dump .= '*';
        }
        $exon_dump .= '  [' . join(' ', @$ex) . "]\n";
      }
      throw("Unexpected: received exon start less than end:\n$exon_dump\n");
    }

    my $exon_len = $exon_end - $exon_start + 1;
    my $cdna_end = $cdna_start + $exon_len - 1;

    if(!defined($transcript_seq_region)) {
      $transcript_seq_region = $seq_region;
    }
    elsif($transcript_seq_region ne $seq_region) {
      push_err(ErrCode::TRANSCRIPT_SCAFFOLD_SPAN);
      ### TBD can probably split transcript rather than discarding
      return ();
    }

    if(!defined($transcript_strand)) {
      $transcript_strand = $exon_strand;
    }
    elsif($transcript_strand != $exon_strand) {
      push_err(ErrCode::TRANSCRIPT_STRAND_FLIP);
      ### TBD can probably split transcript rather than discarding
      return ();
    }

    #
    # calculate the start & end phases
    #

    my ($start_phase, $end_phase);

    if(defined($cur_phase)) {
      if($cdna_start <= $cdna_coding_end) {
        $start_phase = $cur_phase; #start phase is last exons end phase
      } else {
        #the end of the last exon was the end of the CDS
        $start_phase = -1;
      }
    } else {
      #sanity check
      if($cdna_start > $cdna_coding_start && $cdna_start < $cdna_coding_end) {
        throw("Unexpected.  Start of CDS is not in exon?\n" .
              "  exon_cdna_start = $cdna_start\n" .
              "  cdna_coding_start = $cdna_coding_start");
      }
      if($cdna_coding_start == $cdna_coding_start) {
        $start_phase = 0;
      } else {
        $start_phase = -1;
      }
    }

    if($cdna_end < $cdna_coding_start ||
       $cdna_end > $cdna_coding_end) {
      #the end of this exon is outside the CDS
      $end_phase = -1;
    } else {
      #the end of this exon is in the CDS
      #figure out how much coding sequence in the exon
      my $coding_start =
        ($cdna_coding_start > $cdna_start) ? $cdna_coding_start : $cdna_start;
      my $coding_len = $cdna_end - $cdna_start + 1;

      if($start_phase > 0) {
        $coding_len += $start_phase;
      }

      $end_phase = $coding_len % 3;
    }

    my $slice = $slice_adaptor->fetch_by_region('scaffold', $seq_region);

    my $exon = Bio::EnsEMBL::Exon->new
      (-START     => $exon_start,
       -END       => $exon_end,
       -STRAND    => $exon_strand,
       -PHASE     => $start_phase,
       -END_PHASE => $end_phase,
       -SLICE     => $slice);

    $transcript->add_Exon($exon);

    #
    # see if this exon is the start or end exon of the translation
    #

    if($cdna_start <= $cdna_coding_start &&
       $cdna_end   >= $cdna_coding_start) {
      my $translation_start = $cdna_coding_start - $cdna_start + 1;
      $translation->start_Exon($exon);
      $translation->start($translation_start);
    }
    if($cdna_start <= $cdna_coding_end &&
       $cdna_end   >= $cdna_coding_end) {
      my $translation_end = $cdna_coding_end - $cdna_start + 1;
      $translation->end_Exon($exon);
      $translation->end($translation_end);
    }

    $cdna_start = $cdna_end + 1;
    $cur_phase = ($end_phase >= 0) ? $end_phase : undef;
  }

  $transcript->translation($translation);
  push @chimp_transcripts, $transcript;

  return @chimp_transcripts;
}


###############################################################################
# debug
#
###############################################################################

sub debug {
  my $msg  = shift;
  print STDERR "$msg\n" if($verbose);
}




###############################################################################
# usage
#
###############################################################################

sub usage {
  my $msg = shift;

  print STDERR "$msg\n\n" if($msg);

   print STDERR <<EOF;
usage:   perl human2chimp <options>

options: -hdbname <dbname>      human database name

         -hhost <hostname>      human host name (default localhost)

         -huser <user>          human mysql db user with read priveleges

         -hpass <password>      human mysql user password (default none)

         -hport <port>          human mysql db port (default 3306)

         -hassembly <assembly>  human assembly version (default NCBI34)

         -cdbname <dbname>      chimp database name

         -chost <hostname>      chimp host name (default localhost)

         -cuser <user>          chimp mysql db user with read priveleges

         -cpass <password>      chimp mysql user password (default none)

         -cport <port>          chimp mysql db port (default 3306)

         -cassembly <assembly>  chimp assembly version (default BROAD1)

         -help                  display this message

example: perl human2chimp.pl -hdbname homo_sapiens_core_20_34b -hhost ecs2d \\
                             -huser ensro -cdbname pan_trogldytes_core_20_1 \\
                             -chost ecs4 -cport 3350 -cuser ensro

EOF

  exit;
}


