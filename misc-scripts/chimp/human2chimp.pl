use strict;
use warnings;

use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;


my $MAX_CODING_INDEL     = 12; #max allowed coding indel in exon
my $MAX_FRAMESHIFT_INDEL = 5; #max allowed frameshifting indel in exon


my ($hhost, $hdbname, $huser, $hpass, $hport, $hassembly,  #human vars
    $chost, $cdbname, $cuser, $cpass, $cport, $cassembly,  #chimp vars
    $help,  $verbose);

GetOptions('hhost=s'   => \$hhost,
           'hdbname=s' => \$hdbname,
           'huser=s'   => \$huser,
           'hpass=s'   => \$hpass,
           'hport=i'   => \$hport,
           'hassembly=s' => \$hassembly,
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

my $slices = $slice_adaptor->fetch_all('chromosome', $hassembly);

foreach my $slice (@$slices) {

  debug("Chromosome: " . $slice->seq_region_name());

  my $genes = $gene_adaptor->fetch_all_by_Slice($slice);

  foreach my $gene (reverse @$genes) {
    my $transcripts = $gene->get_all_Transcripts();

    foreach my $transcript (@$transcripts) {

      my @transcripts = transfer_transcript($human_db,$chimp_db,$transcript);

    }
  }
}


sub transfer_transcript {
  my $human_db   = shift;
  my $chimp_db   = shift;
  my $transcript = shift;

  my $cs_adaptor      = $human_db->get_CoordSystemAdaptor();
  my $asmap_adaptor   = $human_db->get_AssemblyMapperAdaptor();

  my $chimp_cs = $cs_adaptor->fetch_by_name('scaffold',  $cassembly);
  my $human_cs = $cs_adaptor->fetch_by_name('chromosome', $hassembly);

  my $mapper = $asmap_adaptor->fetch_by_CoordSystems($chimp_cs, $human_cs);

  my $trans_mapper = $transcript->get_TranscriptMapper();

  my $slice_adaptor = $chimp_db->get_SliceAdaptor();

  my $exons = $transcript->get_all_Exons();
  my $translation = $transcript->translation();

  my @chimp_transcripts;

  my $chimp_transcript = Bio::EnsEMBL::Transcript->new();
  my $chimp_translation;

  if($translation) { #watch out for pseudogenes
    $chimp_translation = Bio::EnsEMBL::Translation->new();
    $chimp_translation->start($translation->start());
    $chimp_translation->end($translation->end());
  }

  my $chimp_exon = undef;

 EXON:
  foreach my $exon (@$exons) {

    my @coords = $mapper->map($exon->seq_region_name, $exon->seq_region_start,
                              $exon->seq_region_end, $exon->seq_region_strand,
                              $human_cs);

    if(@coords == 1) {
      my $c = $coords[0];
      if($c->isa('Bio::EnsEMBL::Mapper::Gap')) {
        #complete failure to map exon
        #if this is a UTR exon, this is ok
        my @cds_coords = $trans_mapper->genomic2cds($c->start(),$c->end(),
                                                    $exon->strand);

        #check if this exon was entirely UTR
        if(@cds_coords == 1 &&
           $cds_coords[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
          debug('entire utr exon deletion - ok');

          # we can simply throw out this exon, it was all UTR

          next EXON;
        } else {
          ### TBD: can do more sensible checking and maybe split transcript
          debug('entire coding exon deletion - discarding transcript');
          return @chimp_transcripts;
        }
      } else {
        #exon mapped completely

        debug('completely mapped exon');

        ###TBD maybe should not assume that human phases are set correctly?

        my $slice = $slice_adaptor->fetch_by_region($chimp_cs->name, $c->id,
                                                    undef,undef,undef,
                                                    $chimp_cs->version);


        ### TBD handle problems with seq_region/strand jumping exons...

        $chimp_transcript->add_Exon( Bio::EnsEMBL::Exon->new
                                     ( -START     => $c->start,
                                       -END       => $c->end,
                                       -STRAND    => $c->strand,
                                       -PHASE     => $exon->phase,
                                       -END_PHASE => $exon->end_phase,
                                       -ANALYSIS  => $exon->analysis,
                                       -SLICE     => $slice) );
      }
    } else {
      my $num = scalar(@coords);

      my $exon_position = 0;

      my $chimp_exon = Bio::EnsEMBL::Exon->new();
      my $min_start = undef;
      my $max_end   = undef;
      my $seq_region = undef;

      my $is_start_exon = 0;
      my $is_end_exon   = 0;

      # Fix translation start/end if it needs adjusting here
      if($translation) {

        if($translation->start_Exon()->hashkey() eq $exon->hashkey()) {
          $is_start_exon = 1;
          $chimp_translation->start_Exon($chimp_exon);
        }

        if($translation->end_Exon()->hashkey eq $exon->hashkey()) {
          $is_end_exon = 1;
          $chimp_translation->end_Exon($chimp_exon);
        }
      }


    COORD:
      for(my $i=0; $i < $num; $i++) {
        my $c = $coords[$i];

        $exon_position += $c->length();

        if($c->isa('Bio::EnsEMBL::Mapper::Gap')) {
          #
          # this is a DELETION
          #

          # convert to CDS coords and categorize them
          # gaps are UTR deletions, coords are CDS deletions
          #
          my @cds_coords = $trans_mapper->genomic2cds($c->start(),$c->end(),
                                                      $exon->strand());

          my $cds_delete;
          my @utr_deletes;
          foreach my $cds_coord (@cds_coords) {
            if($cds_coord->isa('Bio::EnsEMBL::Mapper::Gap')) {
              push @utr_deletes, $cds_coord;
            } else {
              #sanity check:
              if($cds_delete) {
                throw("Unexpected: multiple cds deletes for single exon.");
              }
              $cds_delete = $cds_coord;
            }
          }


          if(@utr_deletes > 2) {
            throw("Unexpected: >2 UTR deletes for single exon.");
          }


          foreach my $utr_delete (@utr_deletes) {
            # the UTR deletions may require adjustments to the translation
            debug("utr deletion - ok");

            # Fix translation start/end if it needs adjusting here
            if($is_start_exon) {

              #need to adjust the start of translation if deletion was
              #in 5prime UTR
              if($exon_position < $translation->start()) {
                $chimp_translation->start($chimp_translation->start() -
                                          $utr_delete->length());

                #also need to adjust end translation if this is also the
                #end exon
                if($is_end_exon) {
                  $chimp_translation->end($chimp_translation->end() -
                                          $utr_delete->length());

                }
              }
            }
          }

          if($cds_delete) {
            # if this is a long deletion of coding sequence, then we should
            # throw out this transcript

            my $del_len = $cds_delete->length();
            my $cds_len = get_cds_len($transcript);

            if($cds_len == $del_len) {
              debug("entire cds deletion - discarding transcript");
              return @chimp_transcripts;
            }

            if($del_len > $MAX_CODING_INDEL) {
              debug("long cds deletion ($del_len) - discarding transcript");

              # TBD split transcript in two parts and keep both halves

              return @chimp_transcripts;
            }

            #adjust the translation end accordingly
            if($is_end_exon) {
              #adjust translation end accordingly
              $chimp_translation->end($chimp_translation->end - $del_len);
            }

            my $frame_shift = $del_len %3;

            # if this deletion is short and %3 == 0, then we do not need to
            # do much, but may have to adjust coding end of translation
            if($frame_shift == 0) {
              debug("short ($del_len) in frame cds deletion - ok");

              if($is_end_exon) {
                $chimp_translation->end($chimp_translation->end() - $del_len);
              }

              next COORD;
            }

            if($del_len > $MAX_FRAMESHIFT_INDEL) {
              debug("medium ($del_len) frameshifting cds deletion - " .
                    "discarding transcript");

              #TBD split transcript in two parts and keep both halves

              return @chimp_transcripts;
            }

            # this is a short frameshifting deletion

            # check if deletion is at beginning of CDS
            if($cds_delete->start() == 1) {
              # we need to adjust the start of the exon now to
              # put everything back into frame

              #sanity check:
              if(!$is_start_exon) {
                throw("Unexpected: start of cds is not in start exon");
              }

              my $ex_start = $chimp_exon->start();

              #if there is already UTR, just shift it
              if($chimp_translation->start > 1) {
                debug("frameshift deletion at start of cds - extending utr");
                $chimp_translation->start($chimp_translation->start +
                                          3 - $frame_shift);
              } else {
                $ex_start = $chimp_exon->start();

                #there was no UTR so adjust exon start, leave w/ no UTR
                if(!defined($ex_start)) {
                  #sanity check:
                  if(@coords < $i+2) {
                    throw("Unexpected: no exon start defined and " .
                          "no more coords available");
                  }

                  $ex_start = $coords[$i+1]->start(); #peek at next coord
                }

                debug("frameshift deletion at start of cds - " .
                      "moving exon start");

                $chimp_exon->start($ex_start + 3 - $frame_shift);
              }

              next COORD;
            }

            # check if deletion is at end of CDS
            if($cds_coord->end() == get_cds_len($transcript)) {
              # sanity check:
              if(!$is_end_exon) {
                throw("Unexpected: end of cds is not in end exon");
              }

              # if there was already a UTR just shift it
              if($translation->end() < $exon->length()) {
                $chimp_translation->end($chimp_translation->end() 
                                        - 3 + $frame_shift);
              } else {
                # no utr and at end of CDS, so must be end of exon
                # sanity checks:
                if(!defined($chimp_exon->end())) {
                  throw("Unexpected: at end of exon and no end defined");
                }

                debug("frameshift deletion at end of cds - " .
                      "moving exon end and translation end");
                $chimp_exon->end($chimp_exon->end() - 3 + $frame_shift);
                $chimp_translation->end($chimp_exon->length());
              }

              next COORD;
            }

            # the deletion was not at either end of the cds, but could still
            # be at the beginning or end of the exon

            # check if this is a deletion at the exon start
            if($i == 0) {
              # sanity check:
              if($is_start_exon) {
                throw("Unexpected: delete at start of start exon is " .
                      "not at cds start");
              }

              debug("frameshift deletion at start of exon - " .
                    "adjusting start of exon");
              $chimp_exon->start($chimp_exon->start + 3 - $frameshift);
              if($is_end_exon) {
                debug("adjusting translation end");
                $chimp_exon->translation_end($chimp_exon_translation->end -
                                             3 + $frameshift);
              }

              next COORD;
            }

            #check if this is a deletion at the exon end
            if(@coords == $i+1) {
              # sanity check:
              if($is_end_exon) {
                throw("Unexpected delete at end of end exon is " .
                      "not at cds start");
              }

              debug("frameshift deletion at end of exon - " .
                    "adjusting end of exon");
              $chimp_exon->end($chimp_exon->end - 3 + $frameshift);

              next COORD;
            }


            # now we know that this is a deletion in the middle of the
            # cds

            # for this to be possible there should already be an exon start/end
            # because we should have already seen some exon coords

            # sanity check
            if(!defined($chimp_exon->start())) {
              throw("Unexpected: mid cds delete and no exon start defined");
            }
            if(!defined($chimp_exon->end())) {
              throw("Unexpected: mid cds delete and no exon end defined");
            }

            # split this exon in two parts and make a 'frameshift' intron
            # in the middle

            debug("short ($del_len) frameshifting mid cds deletion - " .
                  "creating fake intron");

            my $ex_len = $chimp_exon->length();

            $chimp_exon->end($chimp_exon->end() - 3 + $frameshift)

            push @chimp_exons, $chimp_exon;

            $chimp_exon = Bio::EnsEMBL::Exon->new
              (-START  => $chimp_exon->end() + 4 - $frame_shift,
               -STRAND => $chimp_exon->strand());

          }

        } else {
          #if first is coord no need to do anything
          next COORD if($i == 0);

          my $prev_c = $coords[$i-1];

          next if($prev_c->isa('Bio::EnsEMBL::Mapper::Gap'));

          #
          # this is an INSERT
          #

          if($c->id() ne $prev_c->id()) {
            debug("exon split across scaffolds - discarding transcript");

            #TBD it may not be necessary to discard transcript if this is UTR

            return @chimp_transcripts;
          }

          if($c->strand() != $prev_c->strand()) {
            debug("exon on both strands - discarding transcript");

            #TBD it may not be necessary to discard transcript if this is UTR

            return @chimp_transcripts;
          }


          # figure out if this insert is in a coding region, and how
          # many bases will be inserted
          # the insert will be in a coding region if it is between 2 coding
          # bases

          my $insert_len;
          my @cds_coords;

          if($c->strand() == 1) {
            my $left_coord  = $exon_position + $exon->start - 1;
            my $right_coord = $left_coord + 1;
            $insert_len = $c->start - $prev_c->end - 1;

            @cds_coords =
              $trans_mapper->genomic2cds($left_coord,$right_coord,1);


          } else {
            my $right_coord = $exon->end() - $exon_position + 1;
            my $left_coord  = $right_coord - 1;
            $insert_len = $prev_c->start() - $c->end() - 1;

            @cds_coords =
              $trans_mapper->genomic2cds($left_coord, $right_coord, -1);
          }

          if($insert_len < 1) {
            debug("inverted exon - discarding transcript");

            #TBD - is there something sensible that can be done here?

            return @chimp_transcripts;
          }


          if(@cds_coords == 1 && 
             $cds_coords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {

            # This is an insert in the CDS

            if($insert_len > $MAX_CODING_INDEL) {
              debug("long cds insert ($insert_len) - discarding transcript");

              # TBD may be ok to keep transcript and split in two?

              return @chimp_transcripts;
            }

            if($insert_len % 3 == 0) {
              #this insert maintains reading frame so ok
              debug("short ($insert_len) in frame cds insert - keeping exon");

              #TBD

              next COORD;
            }

            if($insert_len > $MAX_FRAMESHIFT_INDEL) {
              debug("medium ($insert_len) frameshifting cds insert ".
                    "- discarding transcript");

              #TBD may be ok to keep transcript and split in two

              return @chimp_transcripts;
            }

            debug("short ($insert_len) frameshifting cds insert " .
                  "- creating fake intron");

            # TBD - create fake intron

            next COORD;
          }

          if(@cds_coords == 2) {
            debug("insert b/w cds and utr boundary - ok");
          } elsif(@cds_coords == 1) {
            debug("insert in utr - ok");
          } else {
            throw("Got ", scalar(@cds_coords), " cds_coords, expected 1 or 2");
          }

          #this is an insert in the UTR which is ok

          #TBD may have to adjust translation start/end

        }
      }
    }
  } # foreach exon



  #TBD validation of newly created transcript
  # there are some possible gotchas:
  #  * exons on different scaffolds
  #  * exons on different strands
  #  * exons in which the splice order changed
  #  * transcripts which lost all exons
  # maybe phases should be set here

}



sub get_cds_len {
  my $transcript = shift;

  my $len = 0;

  foreach my $exon (@{$transcript->get_all_translateable_Exons()}) {
    $len += $exon->length();
  }

  return $len;
}

sub debug {
  my $msg  = shift;
  print STDERR "$msg\n" if($verbose);
}


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


