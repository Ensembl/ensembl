# Ensembl module for Bio::EnsEMBL::DBSQL::Utils
#
# Cared for by EnsEMBL (www.ensembl.org)
#
# Copyright GRL and EBI
#
# You may distribute this module under the same terms as perl itself
# POD documentation - main docs before the code



=head1 NAME

Bio::EnsEMBL::DBSQL::Utils - Module having the fset2transcript subroutines

=head1 SYNOPSIS

    use Bio::EnsEMBL::DBSQL::Utils;

    &Bio::EnsEMBL::DBSQL::Utils::fset2transcript($fset_id);

=head1 DESCRIPTION

Module containing the sub routine fset2transcript, 
which creates transcripts from features


=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk


=head1 CONTACT

Ensembl - ensembl-dev@ebi.ac.uk

=head1 APPENDIX

The rest of the documentation details each of the object methods. Internal methods are usually preceded with a _

=cut



package Bio::EnsEMBL::DBSQL::Utils;
use Bio::EnsEMBL::DBSQL::Obj;
use strict;



sub fset2transcript {
    my ($genscan,$contig)=@_;

  
    unless ($genscan->isa ("Bio::EnsEMBL::SeqFeatureI"))
    {print "$genscan must be Bio::EnsEMBL::SeqFeatureI\n";}
    unless ($contig->isa ("Bio::EnsEMBL::DB::ContigI"))
    {print "$contig must be Bio::EnsEMBL::DB::ContigI\n";}

    
    my $transcript = new Bio::EnsEMBL::Transcript;
    $transcript->temporary_id($contig->id . "." . $genscan->raw_seqname);
        
    my @exons;
    my $count= 1;
    
    foreach my $f ($genscan->sub_SeqFeature) {
	
	my $exon  = new Bio::EnsEMBL::Exon;
	$exon->contig_id($contig->internal_id);
	$exon->stable_id($f->id);
	$exon->start    ($f->start);
	$exon->start    ($f->start);
	$exon->end      ($f->end  );
	$exon->strand   ($f->strand);
	$exon->phase    ($f->phase);

	$exon->attach_seq($contig->primary_seq);
	
	push(@exons,$exon);
	$count++;
	
    }
    
    if( $count == 1 ) {
	$genscan->throw("Got a 0 exon genscan");
    }

    my $translation = new Bio::EnsEMBL::Translation;
    #
    # This code got changed due to Translation convention changing. Should work...
    #
    
    if ($exons[0]->strand == 1) {
	@exons = sort {$a->start <=> $b->start} @exons;
    } else {
	@exons = sort {$b->start <=> $a->start} @exons;
    }
    
    $translation->start(1);
    $translation->end($exons[scalar(@exons)-1]->length);
    
    $translation->start_exon($exons[0]);
    $translation->end_exon($exons[$#exons]);
    
    my $endphase = $exons[0]->phase;
    
    foreach my $exon (@exons) {
	
	$exon->phase         ($endphase);
	$transcript->add_Exon($exon);
	$endphase = $exon->end_phase;
	
    }
    
    $transcript->translation($translation);
    
    return $transcript;
}

sub fset2transcript_guess_phases {
    my ($fset,$contig) = @_;

    my $transcript = new Bio::EnsEMBL::Transcript;

    $transcript->temporary_id($contig->id . "." . $fset->id);


    my @exons;
    my $count    = 1;

    foreach my $f ($fset->sub_SeqFeature) {

	my $exon  = new Bio::EnsEMBL::Exon;
	$exon->contig_id($contig->id);
	$exon->start    ($f->start);
	$exon->end      ($f->end  );
	$exon->strand   ($f->strand);
	$exon->attach_seq($contig->primary_seq);
	$exon->phase($f->phase); 
	push(@exons,$exon);
	$count++;
	
    }
	
    my $translation = new Bio::EnsEMBL::Translation;
	
    if ($exons[0]->strand == 1) {
	@exons = sort {$a->start <=> $b->start} @exons;
    } else {
	@exons = sort {$b->start <=> $a->start} @exons;
    }
	
    $translation->start        (1);
    $translation->end          ($exons[$#exons]->end - $exons[$#exons]->start + 1);
    $translation->start_exon($exons[0]);
    $translation->end_exon($exons[$#exons]);
    $transcript->translation($translation);
    
    my $endphase = 0;
    
    foreach my $exon (@exons) {
	
	$exon      ->phase   ($endphase);
	$transcript->add_Exon($exon);

	$endphase = $exon->end_phase;
	
    }

    if ($transcript->translate->seq !~ /\*/) {
	return $transcript;
    }  	

    $endphase = 1;
    
    foreach my $exon (@exons) {
	$exon->phase($endphase);
	$endphase = $exon->end_phase;
    }

    if ($transcript->translate->seq !~ /\*/) {
	return $transcript;
    }  	

    $endphase = 2;
    
    foreach my $exon (@exons) {
	$exon->phase($endphase);
	$endphase = $exon->end_phase;
    }
    
    if ($transcript->translate->seq !~ /\*/) {
	return $transcript;
    }  	
}

sub fset2transcript_3frame {
  my ($fset,$contig) = @_;

  my @f = $fset->sub_SeqFeature;

  if ($f[0]->strand == 1) {
    @f = sort {$a->start <=> $b->start} @f;
  } else {
    @f = sort {$b->start <=> $a->start} @f;
  }

  my @transcripts;

  my $startphase = 0;

  while ($startphase < 3) {
    my $endphase = $startphase;

    my $transcript = new Bio::EnsEMBL::Transcript;

    push(@transcripts,$transcript);

    $transcript->temporary_id($contig->id . "." . $endphase);

    my $count    = 1;
    my @exons;

   
    foreach my $f (@f) {

      my $exon  = new Bio::EnsEMBL::Exon;

      push(@exons,$exon);

      $exon->temporary_id ($contig->id . ".$count");
      $exon->contig_id($contig->id);
      $exon->start    ($f->start);
      $exon->end      ($f->end  );
      $exon->strand   ($f->strand);
      $exon->attach_seq($contig);
      $exon->phase    ($endphase);
      $exon->score    ($f->score);
      $exon->p_value  ($f->p_value);
      $endphase = $exon->end_phase;

      $transcript->add_Exon($exon);
      $count++;

#      print STDERR "Added exon " . $exon->start . " " . $exon->end . " " . $exon->strand . " " . $exon->phase . " " . $exon->end_phase . "\n";
    }
       
    my $translation = new Bio::EnsEMBL::Translation;

    $translation->temporary_id($contig->id . "." . $fset->id);
    $translation->start        (1);
    $translation->end          ($exons[$#exons]->end - $exons[$#exons]->start + 1);
    $translation->start_exon($exons[0]);
    $translation->end_exon  ($exons[$#exons]);
    $transcript->translation($translation);

 #  print STDERR "Phase $startphase " . $transcript->translate->seq . "\n";

    $startphase++;
  }
  return @transcripts;
}


1;

