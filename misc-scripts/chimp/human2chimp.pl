use strict;
use warnings;

use Getopt::Long;

use Bio::EnsEMBL::DBSQL::DBAdaptor;


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

  foreach my $gene (@$genes) {
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

  my @chimp_exons;

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

          # TBD

          next EXON;
        } else {
          ### TBD: can do more sensible checking and maybe split transcript
          debug('entire coding exon deletion - discarding transcript');
          return ();
        }
      } else {
        #exon mapped completely

        debug('completely mapped exon');

        ###TBD maybe should not assume that human phases are set correctly?

        my $slice = $slice_adaptor->fetch_by_region($chimp_cs->name, $c->id,
                                                    undef,undef,undef,
                                                    $chimp_cs->version);

        push @chimp_exons, Bio::EnsEMBL::Exon->new
          (-START     => $c->start,
           -END       => $c->end,
           -PHASE     => $exon->phase,
           -END_PHASE => $exon->end_phase,
           -ANALYSIS  => $exon->analysis,
           -SLICE     => $slice);
      }
    } else {
      my $num = scalar(@coords);

      my $exon_position = 0;

    COORD:
      for(my $i=0; $i < $num; $i++) {
        my $c = $coords[$i];
        $exon_position += $c->length();

        if($c->isa('Bio::EnsEMBL::Mapper::Gap')) {
          #
          # this is a DELETION
          #

          # convert to CDS coords
          my @cds_coords = $trans_mapper->genomic2cds($c->start(),$c->end(),
                                                      $exon->strand());

          @cds_coords = 
            grep {!$_->isa('Bio::EnsEMBL::Mapper::Gap')} @cds_coords;

          if(@cds_coords == 0) {
            #this deletion was entirely in the UTR, accept it

            debug("utr deletion - ok");

            # TBD fix translation start/end if it needs adjusting here

            next COORD;
          }

          if(@cds_coords != 1) {
            # it should not be possible to get more than one cds coord
            # back since we are dealing with a single exon
            throw("Unexpected, got multiple cds coords for single exon");
          }

          my $cds_coord = $cds_coords[0];

          # if this is a long deletion of coding sequence, then we should
          # throw out this transcript

          my $del_len = $cds_coord->length();

          if($del_len > $MAX_CODING_INDEL) {
            debug("long cds deletion ($del_len) - discarding transcript");

            # TBD split transcript in two parts and keep both halves

            return ();
          }

          # if this deletion is short and %3 == 0, then we do not need to
          # do much, but may have to adjust coding start/end of translation
          if($del_len % 3 == 0) {
            debug("short ($del_len) in frame cds deletion - keeping exon");

            #TBD

            next COORD;
          }

          if($del_len > $MAX_FRAMESHIFT_INDEL) {
            debug("medium ($del_len) frameshifting cds deletion - " .
                  "discarding transcript");

            #TBD split transcript in two parts and keep both halves

            return ();
          }

          # this is a short frameshifting deletion
          # create a 'frameshift intron' and move on
          debug("short ($del_len) frameshifting cds deletion - " .
                "creating fake intron");

          # if this exon the first coding exon, then we need to adjust the
          # coding start of the translation, and possibly the exon start
          # as well

          # if this exon is the last coding exon, then we need to adjust the
          # coding end


        } else {
          #if first is coord no need to do anything
          next COORD if($i == 0);

          my $prev_c = $coords[$i];

          next if($prev_c->isa('Bio::EnsEMBL::Mapper::Gap'));

          #
          # this is an INSERT
          #

          if($c->id() ne $prev_c->id()) {
            debug("exon split across scaffolds - discarding transcript");

            #TBD it may not be necessary to discard transcript if this is UTR

            return ();
          }

          if($c->strand() != $prev_c->strand()) {
            debug("exon on both strands - discarding transcript");

            #TBD it may not be necessary to discard transcript if this is UTR

            return ();
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


          if(@cds_coords == 1 && 
             $cds_coords[0]->isa('Bio::EnsEMBL::Mapper::Coordinate')) {

            # This is an insert in the CDS

            if($insert_len > $MAX_CODING_INDEL) {
              debug("long cds insert ($insert_len) - discarding transcript");

              # TBD may be ok to keep transcript and split in two?

              return ();
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

              return ();
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


