use strict;
use warnings;

use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;

use InterimTranscript;
use InterimExon;
use StatMsg;
use Deletion;
use Insertion;
use Transcript;

use Bio::EnsEMBL::Utils::Exception qw(throw info verbose);


{                               #block to avoid namespace pollution

  my ($hhost, $hdbname, $huser, $hpass, $hport, $hassembly, #human vars
      $hchromosome, $hstart, $hend,
      $chost, $cdbname, $cuser, $cpass, $cport, $cassembly, #chimp vars
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

  if ($hchromosome) {
    my $slice = $slice_adaptor->fetch_by_region('chromosome',
                                                $hchromosome,
                                                $hstart, $hend, undef,
                                                $hassembly);
    if (!$slice) {
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
        my $finished_transcripts =
          create_transcripts($interim_transcript, $slice_adaptor);

        foreach my $ftrans (@$finished_transcripts) {
          if($ftrans->translation()) {
            my $pep = $ftrans->translate->seq();
            print STDERR "\n\n$pep\n\n";

            # sanity check, if translation is defined we expect a peptide
            if(!$pep) {
              dump_translation($ftrans->translation());
              throw("Unexpected Translation but no peptide");
            }
          } else {
            print STDERR "NO TRANSLATION LEFT\n";
          }
        }
      }
    }
  }
}

sub dump_translation {
  my $tl = shift;

  print STDERR "TRANSLATION\n";

  if(!$tl) {
    print STDERR "  undef\n";
    return;
  }

  if($tl->start_Exon()) {
    print STDERR "  start exon = ", $tl->start_Exon->stable_id, "\n";
  } else {
    print STDERR "  start exon = undef\n";
  }

  if($tl->end_Exon()) {
    print STDERR "  end exon = ", $tl->end_Exon->stable_id, "\n";
  } else {
    print STDERR "  end exon = undef\n";
  }

  if(defined($tl->start())) {
    print STDERR "  start = ", $tl->start(), "\n";
  } else {
    print STDERR "  start = undef\n";
  }

  if(defined($tl->end())) {
    print STDERR "  end = ", $tl->end(), "\n";
  } else {
    print STDERR "  end = undef\n";
  }


  return;
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

  if (!$transcript->translation()) { # watch out for pseudogenes
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
    info("Exon: " . $human_exon->stable_id() . " chr=" . 
         $human_exon->slice->seq_region_name() . " start=". 
         $human_exon->seq_region_start());
    # info("  cdna_pos = $chimp_cdna_pos\n  cdna_exon_start=$cdna_exon_start");

    my $chimp_exon = InterimExon->new();
    $chimp_exon->stable_id($human_exon->stable_id());
    $chimp_exon->cdna_start($cdna_exon_start);

    my @coords = $mapper->map($human_exon->seq_region_name(),
                              $human_exon->seq_region_start(),
                              $human_exon->seq_region_end(),
                              $human_exon->seq_region_strand(),
                              $human_cs);

    if (@coords == 1) {
      my $c = $coords[0];

      if ($c->isa('Bio::EnsEMBL::Mapper::Gap')) {
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
      }
    } else {
      #
      # Case 3 : Exon mapped partially
      #

      get_coords_extent(\@coords, $chimp_exon);

      if ($chimp_exon->fail()) {
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

        for (my $i=0; $i < $num; $i++) {
          my $c = $coords[$i];

          if ($c->isa('Bio::EnsEMBL::Mapper::Gap')) {

            #
            # deletion in chimp, insert in human
            #
            Deletion::process_delete(\$chimp_cdna_pos, $c->length(),
                                     $chimp_exon, $chimp_transcript);

          } else {
            # can end up with adjacent inserts and deletions so need
            # to take previous coordinate, skipping over gaps
            my $prev_c = undef;

            for (my $j = $i-1; $j >= 0 && !defined($prev_c); $j--) {
              if ($coords[$j]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
                $prev_c = $coords[$j];
              }
            }

            if ($prev_c) {

              my $insert_len;
              if ($chimp_exon->strand() == 1) {
                $insert_len = $c->start() - $prev_c->end() - 1;
              } else {
                $insert_len = $prev_c->start() - $c->end() - 1;
              }

              #sanity check:
              if ($insert_len < 0) {
                throw("Unexpected - negative insert " .
                      "- undetected exon inversion?");
              }

              if ($insert_len > 0) {

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
        }  # foreach coord
      }
    }

    $cdna_exon_start = $chimp_cdna_pos + 1; 
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

    if (!defined($seq_region)) {
      $seq_region = $c->id();
    } elsif ($seq_region ne $c->id()) {
      $chimp_exon->fail(1);
      $stat_code |= StatMsg::SCAFFOLD_SPAN;
      $chimp_exon->add_StatMsg(StatMsg->new($stat_code));
      return;
    }

    if (!defined($strand)) {
      $strand = $c->strand();
    } elsif ($strand != $c->strand()) {
      $chimp_exon->fail(1);
      $stat_code |= StatMsg::STRAND_FLIP;
      $chimp_exon->add_StatMsg(StatMsg->new($stat_code));
      return;
    }

    if (!defined($start)) {
      $start = $c->start if(!defined($start));
    } else {
      if ($strand == 1 && $start > $c->start()) {
        $chimp_exon->fail(1);
        $stat_code |= StatMsg::INVERT;
        $chimp_exon->add_StatMsg(StatMsg->new($stat_code));
        return;
      }
      if ($strand == -1 && $start < $c->start()) {
        $chimp_exon->fail(1);
        $stat_code |= StatMsg::INVERT;
        $chimp_exon->add_StatMsg(StatMsg->new($stat_code));
        return;
      }

      if ($start > $c->start()) {
        $start = $c->start();
      }
    }
	
    if (!defined($end)) {
      $end = $c->end();
    } else {
      if ($strand == 1 && $end > $c->end()) {
        $chimp_exon->fail(1);
        $stat_code |= StatMsg::INVERT;
        $chimp_exon->add_StatMsg(StatMsg->new($stat_code));
        return;
      }
      if ($strand == -1 && $end < $c->end()) {
        $chimp_exon->fail(1);
        $stat_code |= StatMsg::INVERT;
        $chimp_exon->add_StatMsg(StatMsg->new($stat_code));
        return;
      }
      if ($c->end > $end) {
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
  my $itranscript   = shift;    # interim transcript
  my $slice_adaptor = shift;

  # set the phases of the interim exons
  Transcript::set_iexon_phases($itranscript);

  # check the exons and split transcripts where exons are bad
  my $itranscripts = Transcript::check_iexons($itranscript);

  my @finished_transcripts;

  foreach my $itrans (@$itranscripts) {
    # if there are any exons left in this transcript add it to the list
    if (@{$itrans->get_all_Exons()}) {
      push @finished_transcripts, Transcript::make_Transcript($itrans,
                                                              $slice_adaptor);
    } else {
      info("Transcript ". $itrans->stable_id . " has no exons left\n");
    }
  }

  return \@finished_transcripts;
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


