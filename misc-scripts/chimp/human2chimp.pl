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
use Insertion;

use Bio::EnsEMBL::Utils::Exception qw(throw info verbose);


{  #block to avoid namespace pollution

  my ($hhost, $hdbname, $huser, $hpass, $hport, $hassembly,  #human vars
      $hchromosome, $hstart, $hend,
      $chost, $cdbname, $cuser, $cpass, $cport, $cassembly,  #chimp vars
      $help, $verbose);

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

  verbose('INFO') if($verbose); # turn on prints of info statements

  usage() if($help);
  usage("-hdbname option is required") if (!$hdbname);
  usage("-cdbname option is required") if (!$cdbname);

  $hport ||= 3306;
  $cport ||= 3306;

  $hdbname ||= 'localhost';
  $cdbname ||= 'localhost';

  $hassembly ||= 'NCBI34';
  $cassembly ||= 'BROAD1';

  info("Connecting to human database");

  my $human_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (-host    => $hhost,
     -dbname => $hdbname,
     -pass   => $hpass,
     -user   => $huser,
     -port   => $hport);

  info("Connecting to chimp database");

  my $chimp_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (-host    => $chost,
     -dbname  => $cdbname,
     -pass    => $cpass,
     -user    => $cuser,
     -port    => $cport);


  $human_db->dnadb($chimp_db);

  my $slice_adaptor   = $human_db->get_SliceAdaptor();
  my $gene_adaptor    = $human_db->get_GeneAdaptor();


  info("Fetching chromosomes");

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


  my $cs_adaptor      = $human_db->get_CoordSystemAdaptor();
  my $asmap_adaptor   = $human_db->get_AssemblyMapperAdaptor();

  my $chimp_cs = $cs_adaptor->fetch_by_name('scaffold',  $cassembly);
  my $chimp_ctg_cs = $cs_adaptor->fetch_by_name('contig');
  my $human_cs = $cs_adaptor->fetch_by_name('chromosome', $hassembly);

  my $mapper = $asmap_adaptor->fetch_by_CoordSystems($chimp_cs, $human_cs);


  foreach my $slice (@$slices) {

    info("Chromosome: " . $slice->seq_region_name());
    info("Fetching Genes");

    my $genes = $gene_adaptor->fetch_all_by_Slice($slice);

    foreach my $gene (reverse @$genes) {
      info("Gene: ".$gene->stable_id);
      my $transcripts = $gene->get_all_Transcripts();

      foreach my $transcript (@$transcripts) {
        next if(!$transcript->translation); #skip pseudo genes

        my $interim_transcript = transfer_transcript($transcript, $mapper,
						    $human_cs);
	my @transcripts = 
	  create_transcripts($interim_transcript, $slice_adaptor);
      }
    }
  }
}


###############################################################################
#
# transfer_transcript
#
###############################################################################

sub transfer_transcript {
  my $transcript = shift;
  my $mapper = shift;
  my $human_cs = shift;

  info("Transcript: " . $transcript->stable_id());

  my $human_exons = $transcript->get_all_Exons();

  if(!$transcript->translation()) { # watch out for pseudogenes
    info("pseudogene - discarding");
    return;
  }

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
    $chimp_exon->stable_id($human_exon->stable_id());
    $chimp_exon->cdna_start($cdna_exon_start);

    my @coords = $mapper->map($human_exon->seq_region_name(),
			      $human_exon->seq_region_start(),
                              $human_exon->seq_region_end(),
			      $human_exon->seq_region_strand(),
                              $human_cs);

    if(@coords == 1) {
      my $c = $coords[0];

      if($c->isa('Bio::EnsEMBL::Mapper::Gap')) {
        #
        # Case 1: Complete failure to map exon
        #

	my $entire_delete = 1;

	Deletion::process_delete(\$chimp_cdna_pos, $c->length(),
			         $chimp_exon,
				 $chimp_transcript, $entire_delete);

	$chimp_exon->fail(1);
	$chimp_transcript->add_Exon($chimp_exon);

      } else {
        #
        # Case 2: Exon mapped perfectly
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
      #
      # Case 3 : Exon mapped partially
      #

      get_coords_extent(\@coords, $chimp_exon);

      if($chimp_exon->fail()) {
        # Failed to obtain extent of coords due to scaffold spanning
        # strand flipping, or exon inversion.
	# Treat this as if the exon did not map at all.

	my $entire_delete = 1;

	Deletion::process_delete(\$chimp_cdna_pos,
					$human_exon->length(),
					$chimp_exon,
					$chimp_transcript, $entire_delete);

      } else {

	my $num = scalar(@coords);

	for(my $i=0; $i < $num; $i++) {
	  my $c = $coords[$i];

	  if($c->isa('Bio::EnsEMBL::Mapper::Gap')) {

	    #
	    # deletion in chimp, insert in human
	    #
	    Deletion::process_delete(\$chimp_cdna_pos, $c->length(),
				     $chimp_exon, $chimp_transcript);

	  } else {
	    # can end up with adjacent inserts and deletions so need
	    # to take previous coordinate, skipping over gaps
	    my $prev_c = undef;

	    for(my $j = $i-1; $j >= 0 && !defined($prev_c); $j--) {
	      if($coords[$j]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
		$prev_c = $coords[$j];
	      }
	    }

	    if($prev_c) {

	      my $insert_len;
	      if($chimp_exon->strand() == 1) {
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

		Insertion::process_insert(\$chimp_cdna_pos, $insert_len,
					  $chimp_exon, $chimp_transcript);

		$chimp_cdna_pos += $insert_len;
	      }
	    }

	    $chimp_cdna_pos += $c->length();
	  }
	}
      }  # foreach coord

      $cdna_exon_start += $chimp_exon->length();
    }

    $chimp_transcript->add_Exon($chimp_exon);

  } # foreach exon

  return $chimp_transcript;
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
  my $itranscript   = shift; # interim transcript
  my $slice_adaptor = shift;

  my $transcript_strand = undef;
  my $transcript_seq_region = undef;

  my $transcript = Bio::EnsEMBL::Transcript->new();
  my $translation = Bio::EnsEMBL::Translation->new();

  my @chimp_transcripts;

  my $cdna_start = 1;
  my $cur_phase  = undef;

  #
  # If all of the CDS was deleted, then we do not want to create a translation
  #
  my $has_translation = 1;
  if($itranscript->cdna_coding_start == $itranscript->cdna_coding_end + 1) {
    my $stat_msg = StatMsg->new(StatMsg::TRANSCRIPT | StatMsg::NO_CDS_LEFT);
    $itranscript->add_StatMsg($stat_msg);
    $has_translation = 0;
  }

  foreach my $iexon (@{$itranscript->get_all_Exons()}) {

    ### TBD do something with failed exons, maybe split transcript here?
    next if($iexon->fail());

    #sanity check:
    if($iexon->end() < $iexon->start()) {
      throw("Unexpected: exon start less than end.");
    }

    my $cdna_end = $cdna_start + $iexon->length() - 1;

    if(!defined($transcript_seq_region)) {
      $transcript_seq_region = $iexon->seq_region();
    }
    elsif($transcript_seq_region ne $iexon->seq_region()) {
      my $stat_msg = StatMsg->new(StatMsg::TRANSCRIPT | StatMsg::SCAFFOLD_SPAN);
      $itranscript->add_StatMsg($stat_msg);
      ### TBD can probably split transcript rather than discarding
      return;
    }

    if(!defined($transcript_strand)) {
      $transcript_strand = $iexon->strand();
    }
    elsif($transcript_strand != $iexon->strand()) {
      my $stat_msg = StatMsg->new(StatMsg::TRANSCRIPT | StatMsg::STRAND_FLIP);
      $itranscript->add_StatMsg($stat_msg);
      ### TBD can probably split transcript rather than discarding
      return;
    }

    #
    # calculate the start & end phases
    #

    my ($start_phase, $end_phase);

    if($has_translation) {

      if(defined($cur_phase)) {
	if($cdna_start <= $itranscript->cdna_coding_end()) {
	  $start_phase = $cur_phase; #start phase is last exons end phase
	} else {
	  #the end of the last exon was the end of the CDS
	  $start_phase = -1;
	}
      } else {
	#sanity check
	if($cdna_start > $itranscript->cdna_coding_start() &&
	   $cdna_start < $itranscript->cdna_coding_end()) {
	  throw("Unexpected.  Start of CDS is not in exon?\n" .
		"  exon_cdna_start = $cdna_start\n" .
		"  cdna_coding_start = ".$itranscript->cdna_coding_start());
	}
	if($cdna_start == $itranscript->cdna_coding_start()) {
	  $start_phase = 0;
	} else {
	  $start_phase = -1;
	}
      }

      if($cdna_end < $itranscript->cdna_coding_start() ||
	 $cdna_end > $itranscript->cdna_coding_end()) {
	#the end of this exon is outside the CDS
	$end_phase = -1;
      } else {
	#the end of this exon is in the CDS
	#figure out how much coding sequence in the exon
	
	my $coding_start;
	if($itranscript->cdna_coding_start() > $cdna_start) {
	  $coding_start = $itranscript->cdna_coding_start();
	} else {
	  $coding_start = $cdna_start;
	}
	my $coding_len = $cdna_end - $cdna_start + 1;
	
	if($start_phase > 0) {
	  $coding_len += $start_phase;
	}

	$end_phase = $coding_len % 3;
      }
    } else {
      #if there is no translation, all phases should be -1
      $start_phase = -1;
      $end_phase   = -1;
    }

    my $slice =
      $slice_adaptor->fetch_by_region('scaffold', $iexon->seq_region);

    my $exon = Bio::EnsEMBL::Exon->new
      (-START     => $iexon->start(),
       -END       => $iexon->end(),
       -STRAND    => $iexon->strand(),
       -PHASE     => $start_phase,
       -END_PHASE => $end_phase,
       -SLICE     => $slice);

    $transcript->add_Exon($exon);

    #
    # see if this exon is the start or end exon of the translation
    #

    if($has_translation) {
      if($cdna_start <= $itranscript->cdna_coding_start() &&
	 $cdna_end   >= $itranscript->cdna_coding_start()) {
	my $translation_start =
	  $itranscript->cdna_coding_start - $cdna_start + 1;
	$translation->start_Exon($exon);
	$translation->start($translation_start);
      }
      if($cdna_start <= $itranscript->cdna_coding_end() &&
	 $cdna_end   >= $itranscript->cdna_coding_end()) {
	my $translation_end = $itranscript->cdna_coding_end() - $cdna_start + 1;
	$translation->end_Exon($exon);
	$translation->end($translation_end);
      }

      $cdna_start = $cdna_end + 1;
      $cur_phase = ($end_phase >= 0) ? $end_phase : undef;
    }
  }


  if($has_translation) {
    $transcript->translation($translation);
  }
  push @chimp_transcripts, $transcript;

  return @chimp_transcripts;
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


