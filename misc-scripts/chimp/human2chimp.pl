use strict;
use warnings;

use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Analysis;

use InterimTranscript;
use InterimExon;
use StatMsg;
use Deletion;
use Insertion;
use Transcript;
use Gene;
use StatLogger;
use StatMsg;

use Utils qw(print_exon print_coords print_translation);

use Bio::EnsEMBL::Utils::Exception qw(throw info verbose warning);


{                               # block to avoid namespace pollution

  my ($hhost, $hdbname, $huser, $hpass, $hport, $hassembly, # human vars
      $hchromosome, $hstart, $hend,
      $chost, $cdbname, $cuser, $cpass, $cport, $cassembly, # chimp vars
      $dhost, $ddbname, $duser, $dpass, $dport, # destination db
      $help, $verbose, $logfile, $store);

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
             'dhost=s'   => \$dhost,
             'ddbname=s' => \$ddbname,
             'duser=s'   => \$duser,
             'dpass=s'   => \$dpass,
             'dport=i'   => \$dport,
             'store'     => \$store,
             'logfile=s' => \$logfile,
             'help'      => \$help,
             'verbose'   => \$verbose);

  verbose('INFO') if($verbose); # turn on prints of info statements

  usage() if($help);
  usage("-hdbname option is required") if (!$hdbname);
  usage("-cdbname option is required") if (!$cdbname);
  usage("-ddbname option is required when -store is specified")
    if($store && !$ddbname);

  $hport ||= 3306;
  $cport ||= 3306;

  $hdbname ||= 'localhost';
  $cdbname ||= 'localhost';

  $hassembly ||= 'NCBI34';
  $cassembly ||= 'BROAD1';


  info("Connecting to chimp database");

  my $chimp_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (-host    => $chost,
     -dbname  => $cdbname,
     -pass    => $cpass,
     -user    => $cuser,
     -port    => $cport);

  info("Connecting to human database");

  my $human_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (-host    => $hhost,
     -dbname => $hdbname,
     -pass   => $hpass,
     -user   => $huser,
     -port   => $hport);


  my $dest_db;

  if($store) {
    $dest_db = Bio::EnsEMBL::DBSQL::DBAdaptor->new
      (-host   => $dhost,
       -dbname => $ddbname,
       -pass   => $dpass,
       -user   => $duser,
       -dnadb  => $chimp_db,
       -port   => $dport);

    my $analysis = Bio::EnsEMBL::Analysis->new(-logic_name => 'ensembl');
    $dest_db->get_AnalysisAdaptor->store($analysis);
  }

  StatMsg::set_logger(StatLogger->new($logfile));

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


  my $cs_adaptor      = $chimp_db->get_CoordSystemAdaptor();
  my $asmap_adaptor   = $chimp_db->get_AssemblyMapperAdaptor();

  my $chimp_cs = $cs_adaptor->fetch_by_name('chromosome',  $cassembly);
  my $human_cs = $cs_adaptor->fetch_by_name('chromosome', $hassembly);

  my $mapper = $asmap_adaptor->fetch_by_CoordSystems($chimp_cs, $human_cs);
  $mapper->max_pair_count( 6_000_000 );
  $mapper->register_all();

  my $total_transcripts = 0;

  foreach my $slice (@$slices) {

    info("Chromosome: " . $slice->seq_region_name());
    info("Fetching Genes");

    my $genes = $gene_adaptor->fetch_all_by_Slice($slice);

    foreach my $gene (@$genes) {
      info("Gene: ".$gene->stable_id);
      my $transcripts = $gene->get_all_Transcripts();

      my @finished;

      foreach my $transcript (@$transcripts) {
        next if(!$transcript->translation); #skip pseudo genes

        print STDERR ++$total_transcripts, "\n";

        my $interim_transcript = transfer_transcript($transcript, $mapper,
                                                     $human_cs);
        my $finished_transcripts =
          create_transcripts($interim_transcript, $chimp_db->get_SliceAdaptor);


        # set the translation stable identifier on the finished transcripts
        foreach my $tr (@$finished_transcripts) {
          if($tr->translation() && $transcript->translation()) {
            $tr->translation->stable_id($transcript->translation->stable_id);
            $tr->translation->version($transcript->translation->version);
          }
        }

        # This method call is optional, just for statistics gathering
        # and extra sanity checking / debug output:
        gen_transcript_stats($finished_transcripts);

        push @finished, @$finished_transcripts;
      }

      if($store) {
        Gene::store_gene($dest_db, $gene, \@finished);
      }
    }
  }
}



#
# method which does some sanity checking of generated transcripts,
# gathers some stats about the completed transcripts and prints the
# created translations to STDERR
#
sub gen_transcript_stats {
  my $finished_transcripts = shift;

  my $transcript_count = @$finished_transcripts;
  my $translation_count = 0;
  my $stop_codons_count = 0;

  if($transcript_count > 1) {
    StatMsg->new(StatMsg::TRANSCRIPT | StatMsg::SPLIT);
  }
  elsif($transcript_count== 0) {
    StatMsg->new(StatMsg::TRANSCRIPT | StatMsg::NO_SEQUENCE_LEFT);
  }

  foreach my $ftrans (@$finished_transcripts) {
    if($ftrans->translation()) {
      $translation_count++;
      my $pep = $ftrans->translate->seq();

      print STDERR "\n\n$pep\n\n";

      if($pep =~ /\*/) {
        $stop_codons_count++;
      }

      # sanity check, if translation is defined we expect a peptide
      if(!$pep) {
        print_translation($ftrans->translation());
        throw("Unexpected Translation but no peptide");
      }
    } else {
      print STDERR "NO TRANSLATION LEFT\n";
    }
  }

  # If there were stop codons in one of the split transcripts
  # report it. Report it as 'entire' if all split transcripts had
  # stops.
  if($stop_codons_count) {
    my $code = StatMsg::TRANSCRIPT | StatMsg::DOESNT_TRANSLATE;
    if($stop_codons_count == $translation_count) {
      $code |= StatMsg::ENTIRE;
    } else {
      $code |= StatMsg::PARTIAL;
    }
    StatMsg->new($code);
  }

  if(!$translation_count) {
    StatMsg->new(StatMsg::TRANSCRIPT | StatMsg::NO_CDS_LEFT);
  }

  if($translation_count) {
    if($stop_codons_count) {
      if($translation_count > $stop_codons_count) {
        StatMsg->new(StatMsg::TRANSCRIPT | StatMsg::TRANSLATES |
                     StatMsg::PARTIAL);
      }
    } else {
      StatMsg->new(StatMsg::TRANSCRIPT | StatMsg::TRANSLATES |
                   StatMsg::ENTIRE);
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

  my $chimp_cdna_pos = 0;
  my $cdna_exon_start = 1;

  my $chimp_transcript = InterimTranscript->new();
  $chimp_transcript->stable_id($transcript->stable_id());
  $chimp_transcript->version($transcript->version());
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
    $chimp_exon->start_phase($human_exon->phase);
    $chimp_exon->end_phase($human_exon->end_phase());

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

        @coords = @{merge_coords(\@coords)};

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
            if($i > 0) {

              my $prev_c = $coords[$i-1];
              if($prev_c->isa('Bio::EnsEMBL::Mapper::Coordinate')) {

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

  print_coords($coords);

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


###############################################################################
# merge_coords
#
# merges adjacent deletions and insertions
#
###############################################################################

sub merge_coords {
  my $coords = shift;

  info("Coords BEFORE Merging:");
  print_coords($coords);

  my @out = ($coords->[0]);

  # start at 2nd coord because we look at two coords at a time
  for(my $i = 1; $i < @$coords; $i++) {
    my $c1 = $coords->[$i-1];
    my $c2 = $coords->[$i];

    if($c2->isa('Bio::EnsEMBL::Mapper::Gap')) {
      push @out, $c2;
      next;
    }

    # look for adjacent insertion and deletion

    if($c1->isa('Bio::EnsEMBL::Mapper::Coordinate')) {
      push @out, $c2;
      next;
    }

    if($i < 2) {
      # if there was no previous coordinate, then there is no insert as of yet
      push @out, $c2;
      next;
    }

    my $prev_coord = $coords->[$i-2];
    my $insert_len;
    if ($c2->strand() == 1) {
      $insert_len = $c2->start() - $prev_coord->end() - 1;
    } else {
      $insert_len = $prev_coord->start() - $c2->end() - 1;
    }

    #sanity check:
    if ($insert_len < 0) {
      throw("Unexpected - negative insert - undetected exon inversion?");
    }

    if($insert_len == 0) {
      push @out, $c2;
      next;
    }

    # adjacent deletion and insertion was found
    my $del_len = $c1->length();

    info("Merging adjacent insert and delete: I=$insert_len D=$del_len");

    if($del_len > $insert_len) {
      # deletion consumes insertion

      # shrink deletion (grow match) by length of insert
      $c1->end($c1->end - $insert_len);

      if($c2->strand() == 1) {
        $c2->start($c2->start() - $insert_len);
      } else {
        $c2->end($c2->end() + $insert_len);
      }

      push @out, $c2;
    } else {
      #take previous gap off of list and take previous coord
      pop(@out);
      $c1 = $out[$#out];

      if($insert_len - $del_len == 0) {
        # if the insert len is 0 we can merge with previous match
        if($c2->strand() == 1) {
          $c1->end($c2->end());
        } else {
          $c1->start($c2->start());
        }
      } else {
        # grow this match (and shrink insert) by length that was consumed by
        # delete and add it to the list
        if($c2->strand() == 1) {
          $c2->start($c2->start() - $del_len);
        } else {
          $c2->end($c2->end() + $del_len);
        }

        push @out, $c2;
      }
    }
  }

  info("Coords AFTER Merging:");
  print_coords(\@out);

  return \@out;
}

###############################################################################
# create_transcripts
#
###############################################################################

sub create_transcripts {
  my $itranscript   = shift;    # interim transcript
  my $slice_adaptor = shift;

  # set the phases of the interim exons
  # Transcript::set_iexon_phases($itranscript);

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

         -store                 flag indicating genes are to be stored in a
                                destination database

         -ddbname <dbname>      destination database name

         -dhost <hostname>      destination host name (default localhost)

         -duser <user>          destination mysql db user with write priveleges

         -dpass <password>      destination mysql user password (default none)

         -dport <port>          desitnation mysql db port (default 3306)

         -verbose               flag indicating debug and progress messages
                                should be printed to STDERR

         -logfile <file>        File to which status messages should be logged
                                to.  The contents of this file can be
                                interpreted using the get_stats.pl script.

         -help                  display this message

example: perl human2chimp.pl -hdbname homo_sapiens_core_20_34b -hhost ecs2d \\
                             -huser ensro -cdbname pan_troglodytes_core_20_1 \\
                             -chost ecs4 -cport 3350 -cuser ensro \\
                             -store -ddbname pt_genes -dhost ecs1d \\
                             -duser ensadmin -dpass secret -logfile out.txt

EOF

  exit;
}




