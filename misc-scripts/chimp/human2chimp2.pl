use strict;
use warnings;

use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Exon;

use Bio::EnsEMBL::Utils::Exception qw(throw);

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


debug("Fetching chromosomes");

my $slices = $slice_adaptor->fetch_all('chromosome', $hassembly);

foreach my $slice (@$slices) {

  debug("Chromosome: " . $slice->seq_region_name());
  debug("Fetching Genes");

  my $genes = $gene_adaptor->fetch_all_by_Slice($slice);

  foreach my $gene (reverse @$genes) {
    debug("Gene: ".$gene->stable_id);
    my $transcripts = $gene->get_all_Transcripts();

    foreach my $transcript (@$transcripts) {

      my @transcripts = transfer_transcript($human_db,$chimp_db,$transcript);

    }
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

  my $cs_adaptor      = $human_db->get_CoordSystemAdaptor();
  my $asmap_adaptor   = $human_db->get_AssemblyMapperAdaptor();

  my $chimp_cs = $cs_adaptor->fetch_by_name('scaffold',  $cassembly);
  my $human_cs = $cs_adaptor->fetch_by_name('chromosome', $hassembly);

  my $mapper = $asmap_adaptor->fetch_by_CoordSystems($chimp_cs, $human_cs);
  my $exons = $transcript->get_all_Exons();

  if(!$transcript->translation()) { #watch out for pseudogenes
    debug("pseudogene - discarding");
    return ();
  }

  my $trans_mapper = $transcript->get_TranscriptMapper();

  my $chimp_cdna_pos = 0;
  my $cdna_exon_start = 1;
  my $cdna_coding_start = $transcript->cdna_coding_start();
  my $cdna_coding_end   = $transcript->cdna_coding_end();

  my @chimp_exons;


 EXON:
  foreach my $exon (@$exons) {

    my @coords = $mapper->map($exon->seq_region_name, $exon->seq_region_start,
                              $exon->seq_region_end, $exon->seq_region_strand,
                              $human_cs);

    if(@coords == 1) {
      my $c = $coords[0];

      if($c->isa('Bio::EnsEMBL::Mapper::Gap')) {
        #
        # complete failure to map exon
        #

        #if this is a UTR exon, this is ok
        my @cds_coords = $trans_mapper->genomic2cds($c->start(),$c->end(),
                                                    $exon->strand);

        #check if this exon was entirely UTR
        if(@cds_coords == 1 &&
           $cds_coords[0]->isa('Bio::EnsEMBL::Mapper::Gap')) {
          debug('entire utr exon deletion - ok');

          # we can simply throw out this exon, it was all UTR

        } else {
          ### TBD: can do more sensible checking and maybe split transcript
          debug('entire coding exon deletion - discarding transcript');
          return ();
        }
      } else {
        #
        # exon mapped completely
        #

        debug('completely mapped exon');

        ### TBD handle problems with seq_region/strand jumping exons...

        $chimp_cdna_pos += $c->length();
        $cdna_exon_start += $c->length();

        push @chimp_exons, [$c->start(), $c->end(), $c->strand(), $c->id()];
      }
    } else {
      my $num = scalar(@coords);

      my $result = get_coords_extent(\@coords);

      if(!$result) {
        #failed to obtain extent of coords due to scaffold spanning
        #strand flipping, or exon inversion
	
        ### TBD - fix this, may be ok to drop exon and continue esp. if exon
        ###       is entirely UTR
        return ();
      }

      my($exon_start, $exon_end, $exon_strand, $seq_region) = @$result;

      for(my $i=0; $i < $num; $i++) {
        my $c = $coords[$i];

        if($c->isa('Bio::EnsEMBL::Mapper::Gap')) {

          #
          # deletion in chimp, insert in human
          #
          my $result =
            process_deletion(\$chimp_cdna_pos, $c->length(),
                             \$exon_start, \$exon_end, $exon_strand,
                             \$cdna_exon_start,
                             \$cdna_coding_start, \$cdna_coding_end);

          return () if(!defined($result));

          push @chimp_exons, @$result;

        } else {
          my $prev_c = $coords[$i-1];

          if($i > 0 && !$prev_c->isa('Bio::EnsEMBL::Mapper::Gap')) {

            #
            # insert in chimp, deletion in human
            #

            my $insert_len;
            if($exon_strand == 1) {
              $insert_len = $c->start() - $prev_c->end() - 1;
            } else {
              $insert_len = $prev_c->start() - $c->end() - 1;
            }

            my $result = process_insertion(\$chimp_cdna_pos, $insert_len,
                                       \$exon_start, \$exon_end, $exon_strand,
                                       \$cdna_exon_start,
                                       \$cdna_coding_start, \$cdna_coding_end);

            return () if(!defined($result));

            push @chimp_exons, @$result;

            $chimp_cdna_pos += $insert_len;

          }
          $chimp_cdna_pos += $c->length()
        }
      }  # foreach coord

      $cdna_exon_start += $exon_end - $exon_start + 1;
    }
  } # foreach exon
}

###############################################################################
# process_deletion
#
# Processes a deletion in an exon by adjusting the exon start/end
# and adjusting the translation start/end of the transcript as necessary
#
# If the deletion is in the CDS and considered to be too long to accomodate
# this method will return undef.
#
# On success this this method will return a listref of extra exons that
# needed to be created as sometimes it is necessary to split exons into
# multiple parts.
#
###############################################################################

sub process_deletion {
  my $chimp_cdna_pos_ref = shift;
  my $del_len            = shift;

  my $exon_start_ref     = shift;
  my $exon_end_ref       = shift;
  my $exon_strand        = shift;

  my $cdna_exon_start_ref = shift;

  my $cdna_coding_start_ref = shift;
  my $cdna_coding_end_ref   = shift;

  my $del_start = $$chimp_cdna_pos_ref + 1;
  my $del_end   = $$chimp_cdna_pos_ref + $del_len;

  my $exon_len = $$exon_end_ref - $$exon_start_ref + 1;
  my $cdna_exon_end = $$cdna_exon_start_ref + $exon_len - 1;

  print STDERR "del_start = $del_start\n";
  print STDERR "del_end   = $del_end\n";

  # sanity check, deletion should be completely in exon boundaries
  if($del_start < $$cdna_exon_start_ref ||
     $del_end   > $cdna_exon_end) {

    # it is ok if this deletion is immediately in front of the exon or
    # immediately behind the exon
    if($del_start != $cdna_exon_end + 1 &&
       $del_end   != $$cdna_exon_start_ref - 1) {

      throw("Unexpected: deletion is outside of exon boundary\n" .
            "     del_start       = $del_start\n" .
            "     del_end         = $del_end\n" .
            "     cdna_exon_start = $$cdna_exon_start_ref\n" .
            "     cdna_exon_end   = $cdna_exon_end");
    }
  }


  #
  # case 1: delete covers entire CDS
  #
  if($del_start <= $$cdna_coding_start_ref &&
     $del_end >= $$cdna_coding_end_ref) {

    # nothing can be done with this exon since entire CDS is deleted
    debug("entire cds deleted - discarding transcript");
    return undef;
  }

  #
  # case 2: delete overlaps UTR and start of CDS
  #
  elsif($del_start < $$cdna_coding_start_ref &&
        $del_end   >= $$cdna_coding_start_ref) {
    my $utr_del_len = $$cdna_coding_start_ref - $del_start;
    my $cds_del_len = $del_len - $utr_del_len;

    debug("delete in 5' utr ($utr_del_len) and start of cds ($cds_del_len)");

    if($cds_del_len > $MAX_CODING_INDEL) {
      # too much coding sequence has been deleted
      debug("cds deletion is too long ($cds_del_len) - " .
            "discarding transcript");
      return undef;
    }

    # move up CDS to account for delete in 5prime UTR
    $$cdna_coding_start_ref -= $utr_del_len;
    $$cdna_coding_end_ref   -= $utr_del_len;

    #move up CDS end to account for delete in CDS
    $$cdna_coding_end_ref   -= $cds_del_len;

    my $frameshift = $cds_del_len % 3;

    if($frameshift) {
      if($cds_del_len > $MAX_FRAMESHIFT_INDEL) {
        # frameshift deletion is too long
        debug("frameshift deletion is too long ($cds_del_len) " .
              "- discarding transcript");
        return undef;
      }

      # move down CDS start to put reading frame back (shrink CDS)
      debug("shifting cds start to restore reading frame");
      $$cdna_coding_start_ref += 3 - $frameshift;

      if($$cdna_coding_start_ref >= $$cdna_coding_end_ref) {
        debug("no cds left after shift - discarding transcript");
        return undef;
      }
    }
  }

  #
  # case 3: delete overlaps end of CDS and UTR
  #
  elsif($del_start <= $$cdna_coding_end_ref &&
        $del_end   > $$cdna_coding_end_ref) {
    my $cds_del_len = $$cdna_coding_end_ref - $del_start + 1;
    my $utr_del_len = $del_len - $cds_del_len;

    debug("delete in 3' utr ($utr_del_len) and end of cds ($cds_del_len)");

    if($cds_del_len > $MAX_CODING_INDEL) {
      # too much coding sequence has been deleted
      debug("cds deletion is too long ($cds_del_len) - " .
            "discarding transcript");
      return undef;
    }

    # move up CDS end to account for CDS deletion
    $$cdna_coding_end_ref -= $cds_del_len;

    my $frameshift = $cds_del_len % 3;

    if($frameshift) {
      if($cds_del_len > $MAX_FRAMESHIFT_INDEL) {
        debug("frameshift deletion is too long ($cds_del_len) " .
              "- discarding transcript");
      return undef;
      }

      #move up CDS end to put reading frame back (shrink CDS)
      debug("shifting cds end to restore reading frame");
      $$cdna_coding_end_ref -= 3 - $frameshift;

      if($$cdna_coding_start_ref >= $$cdna_coding_end_ref) {
        debug("no cds left after shift - discarding transcript");
        return undef;
      }
    }
  }

  #
  # case 4: delete is at start of CDS
  #
  elsif($del_start == $$cdna_coding_start_ref) {
    debug("delete ($del_len) at start of cds");

    if($del_len > $MAX_CODING_INDEL) {
      # too much coding sequence has been deleted
      debug("cds deletion is too long ($del_len) - " .
            "discarding transcript");
      return undef;
    }

    # move up CDS end to account for CDS deletion
    $$cdna_coding_end_ref -= $del_len;

    my $frameshift = $del_len % 3;

    if($frameshift) {
      if($del_len > $MAX_FRAMESHIFT_INDEL) {
        debug("frameshift deletion is too long ($del_len) " .
              "- discarding transcript");
        return undef;
      }

      # move down CDS start to put reading frame back (shrink CDS)
      debug("shifting cds start to restore reading frame");
      $$cdna_coding_start_ref += 3 - $frameshift;

      if($$cdna_coding_start_ref >= $$cdna_coding_end_ref) {
        debug("no cds left after shift - discarding transcript");
        return undef;
      }
    }
  }

  #
  # case 5: delete is at end of CDS
  #
  elsif($del_end == $$cdna_coding_end_ref) {
    debug("delete ($del_len) at end of cds");

    if($del_len > $MAX_CODING_INDEL) {
      # too much coding sequence has been deleted
      debug("cds deletion is too long ($del_len) - " .
            "discarding transcript");
      return undef;
    }

    # move up CDS end to account for CDS deletion
    $$cdna_coding_end_ref -= $del_len;

    my $frameshift = $del_len % 3;

    if($frameshift) {
      if($del_len > $MAX_FRAMESHIFT_INDEL) {
        debug("frameshift deletion is too long ($del_len) " .
              "- discarding transcript");
        return undef;
      }

      # move up CDS end to put reading frame back (shrink CDS)
      $$cdna_coding_end_ref -= 3 - $frameshift;

      if($$cdna_coding_start_ref >= $$cdna_coding_end_ref) {
        debug("no cds left after shift - discarding transcript");
        return undef;
      }
    }
  }

  #
  # case 6: delete is in middle of CDS
  #
  elsif($del_end   > $$cdna_coding_start_ref &&
        $del_start < $$cdna_coding_end_ref) {
    debug("delete ($del_len) in middle of cds");

    if($del_len > $MAX_CODING_INDEL) {
      # too much coding sequence has been deleted
      debug("cds deletion is too long ($del_len) - " .
            "discarding transcript");
      return undef;
    }

    #move up CDS end to account for CDS deletion
    $$cdna_coding_start_ref -= $del_len;

    my $frameshift = $del_len % 3;

    if($frameshift) {
      if($del_len > $MAX_FRAMESHIFT_INDEL) {
        debug("frameshift deletion is too long ($del_len) " .
              " - discarding transcript");
        return undef;
      }

      debug("introducing frameshift intron to maintain reading frame");

      # this is going to require splitting the exon
      # to make a frameshift deletion

      #first exon is going to end right before deletion
      my $first_len  = $del_start - $$cdna_coding_start_ref;
      my $intron_len = 3 - $frameshift;

      #second exon is going to start right after 'frameshift intron'
      $$cdna_exon_start_ref += $first_len + $intron_len;

      #reduce the length of the CDS by the length of the new intron
      $$cdna_coding_end_ref -= $intron_len;

      ### TBD may have to check we have not run up to end of CDS here

      if($exon_strand == 1) {
        #end the current exon at the beginning of the deletion
        my $out_exon = [$$exon_start_ref, $$exon_start_ref + $first_len - 1];
        $$exon_start_ref += $first_len + $intron_len;
        return [$out_exon];

      } else {
        my $out_exon = [$$exon_end_ref - $first_len + 1, $$exon_end_ref];
        $$exon_end_ref   -= $first_len + $intron_len;
        return [$out_exon];
      }
    }
  }

  #
  # case 7: delete is in 5prime UTR
  #
  elsif($del_end < $$cdna_coding_start_ref) {
    debug("delete ($del_len) in 5' utr");

    # just move up the CDS
    $$cdna_coding_start_ref -= $del_len;
    $$cdna_coding_end_ref   -= $del_len;
  }

  #
  # case 8: delete is in 3prime UTR
  #
  elsif($del_start > $$cdna_coding_end_ref) {
    debug("delete ($del_len) in 3' utr");

    #do not have to do anything
  }

  #
  # default : sanity check
  #
  else {
    throw("Unexpected delete case encountered");
  }

  return [];
}

###############################################################################
# process_insertion
#
###############################################################################

sub process_insertion {
  my $chimp_cdna_pos_ref = shift;   #basepair to left of insert

  my $insert_len         = shift;

  my $exon_start_ref     = shift;
  my $exon_end_ref       = shift;
  my $exon_strand        = shift;

  my $cdna_exon_start_ref = shift;

  my $cdna_coding_start_ref = shift;
  my $cdna_coding_end_ref   = shift;

  my $exon_len = $$exon_end_ref - $$exon_start_ref + 1;
  my $cdna_exon_end = $$cdna_exon_start_ref + $exon_len - 1;

  # sanity check, insert should be completely in exon boundaries
  if($$chimp_cdna_pos_ref < $$cdna_exon_start_ref ||
     $$chimp_cdna_pos_ref >= $cdna_exon_end) {
    throw("Unexpected: deletion is outside of exon boundary");
  }


  #
  # case 1: insert in CDS
  #
  if($$chimp_cdna_pos_ref >= $$cdna_coding_start_ref &&
     $$chimp_cdna_pos_ref < $$cdna_coding_end_ref) {

    debug("insertion in cds ($insert_len)");

    if($insert_len > $MAX_CODING_INDEL) {
      # too much coding sequence has been inserted
      debug("cds insertion is too long ($insert_len) - " .
            "discarding transcript");
      return undef;
    }

    #adjust CDS end accordingly
    $$cdna_coding_end_ref += $insert_len;

    my $frameshift = $insert_len % 3;

    if($frameshift) {
      if($insert_len > $MAX_FRAMESHIFT_INDEL) {
        debug("frameshift insertion is too long ($insert_len) - " .
              "discarding transcript");
        return undef;
      }

      # need to create frameshift intron to get reading frame back on track
      # exon needs to be split into two

      debug("introducing frameshift intron to maintain reading frame");

      #first exon ends right before insert
      my $first_len  = $$chimp_cdna_pos_ref - $$cdna_coding_start_ref + 1;

      # frame shift intron eats into start of inserted region
      # second exon is going to start right after 'frameshift intron'
      $$cdna_exon_start_ref += $first_len + $frameshift;

      #reduce the length of the CDS by the length of the new intron
      $$cdna_coding_end_ref -= $frameshift;

      ### TBD may have to check we have not run up to end of CDS here

      if($exon_strand == 1) {
        #end the current exon at the beginning of the deletion
        my $out_exon = [$$exon_start_ref, $$exon_start_ref + $first_len - 1];
        $$exon_start_ref += $first_len + $frameshift;
        return [$out_exon];

      } else {
        my $out_exon = [$$exon_end_ref - $first_len + 1, $$exon_end_ref];
        $$exon_end_ref   -= $first_len + $frameshift;
        return [$out_exon];
      }
    }
  }

  #
  # case 2: insert in 5 prime UTR (or between 5prime UTR and CDS)
  #
  elsif($$chimp_cdna_pos_ref < $$cdna_coding_start_ref) {
    debug("insertion ($insert_len) in 5' utr");

    #shift the coding region down as result of insert
    $$cdna_coding_start_ref += $insert_len;
    $$cdna_coding_end_ref   += $insert_len;
  }

  #
  # case 3: insert in 3 prime UTR (or between 3prime UTR and CDS)
  #
  elsif($$chimp_cdna_pos_ref >= $$cdna_coding_end_ref) {
    debug("insert ($insert_len) in 3' utr");

    #do not have to do anything
  }

  #
  # default: sanity check
  #
  else {
    throw("Unexpected insert case encountered");
  }
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

  my($start, $end, $strand, $seq_region);

  foreach my $c (@$coords) {
    next if($c->isa('Bio::EnsEMBL::Mapper::Gap'));

    if(!defined($seq_region)) {
      $seq_region = $c->id();
    }
    elsif($seq_region ne $c->id()) {
      debug("coords spans multiple scaffolds - unable to get extent");
      return undef;
    }

    if(!defined($strand)) {
      $strand = $c->strand();
    }
    elsif($strand != $c->strand()) {
      debug("coords flip strands - unable to get extent");

      return undef;
    }

    if(!defined($start)) {
      $start = $c->start if(!defined($start));
    } elsif($start > $c->start()) {
      if($strand == 1) {
	debug("coord inversion - unable to get extent");

	return undef;
      } else {
	$strand = $c->strand();
      }
    }
	
    if(!defined($end) || $c->end > $end) {
      $end   = $c->end();
    }
  }

  return [$start, $end, $strand, $seq_region];
}


###############################################################################
# get_cds_len
#
# returns the CDS length of the transcript
#
###############################################################################
sub get_cds_len {
  my $transcript = shift;

  my $len = 0;

  foreach my $exon (@{$transcript->get_all_translateable_Exons()}) {
    $len += $exon->length();
  }

  return $len;
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


