use strict;
use warnings;

package Utils;

use Exporter;
use Bio::EnsEMBL::Utils::Exception qw(info);

use vars qw(@ISA @EXPORT_OK);

@ISA = qw(Exporter);

@EXPORT_OK = qw(print_exon print_coords print_translation print_three_phase_translation);



#
# prints debugging info
#
sub print_exon {
  my $exon = shift;
  my $tr = shift;

  if (!$exon) {
    throw("Exon undefined");
  }

  info(" ".$exon->stable_id());

  info("  cdna_start = ".$exon->cdna_start())
    if(defined($exon->cdna_start()));

  info("  cdna_end   = ". $exon->cdna_end())
    if(defined($exon->cdna_end()));

  info("  start             = ". $exon->start())
    if(defined($exon->start()));

  info("  end               = ". $exon->end())
    if(defined($exon->end()));

  info("  strand            = ". $exon->strand())
    if(defined($exon->strand()));

  if($exon->fail) {
    info("  FAILED");
  }

  if($tr) {
    info(" TRANSCRIPT:");
    info("  cdna_coding_start = ". $tr->cdna_coding_start());
    info("  cdna_coding_end   = ". $tr->cdna_coding_end(). "\n");
  }

  return;
}




sub print_coords {
  my $cs = shift;

  foreach my $c (@$cs) {
    if($c->isa('Bio::EnsEMBL::Mapper::Gap')) {
      info("  GAP ". $c->length());
    } else {
      info("  ". $c->start. '-'. $c->end. ' ('.$c->strand.")");
    }
  }
}


sub print_translation {
  my $tl = shift;

  info("TRANSLATION");

  if(!$tl) {
    info("  undef");
    return;
  }

  if($tl->start_Exon()) {
    info("  start exon = ", $tl->start_Exon->stable_id());
  } else {
    info("  start exon = undef");
  }

  if($tl->end_Exon()) {
    info("  end exon = ", $tl->end_Exon->stable_id);
  } else {
    info("  end exon = undef");
  }

  if(defined($tl->start())) {
    info("  start = ", $tl->start());
  } else {
    info("  start = undef");
  }

  if(defined($tl->end())) {
    info("  end = ", $tl->end());
  } else {
    info("  end = undef");
  }

  return;
}


sub print_three_phase_translation {
  my $transcript = shift;

  return if(!$transcript->translation());

  my $orig_phase = $transcript->start_Exon->phase();

  foreach my $phase (0,1,2) {
    info("======== Phase $phase translation: ");
    $transcript->start_Exon->phase($phase);
    info("Peptide: " . $transcript->translate->seq() . "\n\n===============");
  }

  $transcript->start_Exon->phase($orig_phase);

  return;
}



1;
