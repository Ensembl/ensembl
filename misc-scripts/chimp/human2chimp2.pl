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

      my $result = get_coords_extent(@coords);

      if(!$result) {
	#failed to obtain extent of coords due to scaffold spanning
	#strand flipping, or exon inversion
	
	### TBD - fix this, may be ok to drop exon and continue esp. if exon
	###       is entirely UTR
	return @chimp_transcripts;
      }

      my($exon_start, $exon_end, $exon_strand, $seq_region) = @$result;

      my $human_cdna_pos = 0;
      
      my $cdna_coding_start = $transcript->cdna_coding_start();
      my $cdna_coding_end   = $transcript->cdna_coding_end();

      my $chimp_cdna_pos = 0;

    COORD:
      for(my $i=0; $i < $num; $i++) {
        my $c = $coords[$i];

        if($c->isa('Bio::EnsEMBL::Mapper::Gap')) {
	  #deletion in chimp, insert in human
	  process_deletion(\$human_cdna_pos, \$chimp_cdna_pos, $c->length(),
			   \$exon_start, \$exon_end, $exon_strand,
			   \$cdna_coding_start, \$cdna_coding_end);

}


sub process_deletion {
  my $human_cdna_pos_ref = shift;
  my $chimp_cdna_pos_ref = shift;
  my $del_length         = shift;

  my $exon_start_ref     = shift;
  my $exon_end_ref       = shift;
  my $exon_strand        = shift;

  my $cdna_coding_start_ref = shift;
  my $cdna_coding_end_ref   = shift;

  my $$human_cdna_pos_ref += $del_length;

  my $del_start = $$chimp_cdna_pos_ref + 1;
  my $del_end   = $$chimp_cdna_pos_ref + $del_length;

  # case 1: delete covers entire CDS
  if($del_start <= $$cdna_coding_start_ref &&
     $del_end >= $$cdna_coding_end_ref) {
  }

  # case 2: delete overlaps UTR and start of CDS
  elsif($del_start < $$cdna_coding_start_ref &&
	$del_end   >= $$cdna_coding_start_ref) {
  }

  # case 3: delete overlaps end of CDS and UTR
  elsif($del_start <= $$cdna_coding_end_ref &&
	$del_end   > $$cdna_coding_end_ref) {
  }

  # case 4: delete is at start of CDS
  elsif($del_start == $$cdna_coding_start_ref) {
  }

  # case 5: delete is at end of CDS
  elsif($del_end == $$cdna_coding_end_ref) {
  }

  # case 6: delete is in middle of CDS
  elsif($del_end   > $$cdna_coding_start_ref &&
	$del_start < $$cdna_coding_end_ref) {

  }

  # case 7: delete is in UTR
  elsif($del_end < $$cdna_coding_start_ref ||
	$del_start > $$cdna_coding_end_ref) {

  }

  # default : sanity check
  else {
    throw("Unexpected delete case encountered");
  }

  if($$chimp_cdna_pos_ref < $$cdna_coding_start_ref) {
    $$cdna_coding_start_ref -= $del_length;
    $$cdna_coding_start_ref = 1 if($$cdna_coding_start_ref < 1);
  }

  if($$chimp_cdna_pos_ref < $$cdna_coding_end_ref) {
    $$cdna_coding_end_ref -= $del_length;
    $$cdna_coding_start_ref = 1 if($$cdna_coding_start_ref < 0);
  }


}


sub process_insert {
}


#
# given a list of coords returns the start, end, strand, seq_region
# of the span of the coords
#
# undef is returned if the coords flip strands, have an inversion,
# or cross multiple seq_regions
#
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



#
# returns the CDS length of the transcript
#
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


