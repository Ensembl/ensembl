=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut

=head1 NAME

Bio::EnsEMBL::TranscriptFactory - Module having the fset2transcript*
subroutines

=head1 SYNOPSIS

  use Bio::EnsEMBL::TranscriptFactory;

  &Bio::EnsEMBL::TranscriptFactory::fset2transcript($fset_id);

=head1 DESCRIPTION

Module containing the subroutines fset2transcript*, 
which create transcripts from features (formally housed in
Bio::EnsEMBL::DBSQL::Utils).

=head1 METHODS

=cut

package Bio::EnsEMBL::TranscriptFactory;

use strict;

use Bio::EnsEMBL::Exon;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Transcript;

sub fset2transcript {
    my ($genscan,$contig)=@_;

  
    unless ($genscan->isa ("Bio::EnsEMBL::SeqFeatureI"))
    {print "$genscan must be Bio::EnsEMBL::SeqFeatureI\n";}
     
    my $transcript = new Bio::EnsEMBL::Transcript;
    $transcript->temporary_id($contig->id . "." . $genscan->seqname);
        
    my @exons;
    my $count= 1;
    
    foreach my $f ($genscan->sub_SeqFeature) {
  
	my $exon  = new Bio::EnsEMBL::Exon;
	$transcript->add_Exon($exon);
        $exon->contig   ($contig);
	$exon->start    ($f->start);
	$exon->end      ($f->end  );
	$exon->strand   ($f->strand);
	$exon->phase    ($f->phase);
	$exon->end_phase( ($exon->phase + $exon->length)%3 );
	#$exon->score($f->score);
#	$exon->p_value($f->p_value);
	$exon->slice($contig->primary_seq);
	
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
    
    $translation->start_Exon($exons[0]);
    $translation->end_Exon($exons[$#exons]);

    my $endphase = $exons[0]->end_phase;
    
    foreach my $exon (@exons) {
	
      if ( $exon == $exons[0] ){
	next;
      }
      $exon->phase         ($endphase);
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
        $exon->contig   ($contig);
	$exon->start    ($f->start);
	$exon->end      ($f->end  );
	$exon->strand   ($f->strand);
	#$exon->score($f->score);
#	$exon->p_value($f->p_value);
	$exon->slice($contig);
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
    $translation->start_Exon($exons[0]);
    $translation->end_Exon($exons[$#exons]);
    $transcript->translation($translation);
    
    my $endphase = 0;
    
    foreach my $exon (@exons) {
	
	$exon      ->phase   ($endphase);
	$transcript->add_Exon($exon);

	$endphase = $exon->end_phase(($exon->phase + $exon->length)%3);
	
    }


    if ($transcript->translate->seq !~ /\*/) {
	return $transcript;
    }  	

    $endphase = 1;
    
    foreach my $exon (@exons) {
	$exon->phase($endphase);
	$endphase = $exon->end_phase(($exon->phase + $exon->length)%3);
    }

    if ($transcript->translate->seq !~ /\*/) {
	return $transcript;
    }  	

    $endphase = 2;
    
    foreach my $exon (@exons) {
	$exon->phase($endphase);
	$endphase = $exon->end_phase(($exon->phase + $exon->length)%3);
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
      #print "exon seqname = ".$f->seqname."\n";
      my $exon  = new Bio::EnsEMBL::Exon;
      #print STDERR "exon ".$f->gffstring."\n";
      push(@exons,$exon);
      $exon->seqname($f->seqname);
      $exon->temporary_id ($contig->id . ".$count");
      $exon->contig   ($contig);
      $exon->start    ($f->start);
      $exon->end      ($f->end  );
      $exon->strand   ($f->strand);
      $exon->slice($contig);
      $exon->phase    ($endphase);
      $exon->end_phase( ($exon->phase + $exon->length)%3 );
      #$exon->score    ($f->score);
#      $exon->p_value  ($f->p_value);
      $endphase = $exon->end_phase;

      $transcript->add_Exon($exon);
      $count++;

      #print STDERR "Added exon start " . $exon->start . " end " . $exon->end . " strand " . $exon->strand . " score " . $exon->score . " pvalue " . $exon->p_value . "\n";
    }
       
    my $translation = new Bio::EnsEMBL::Translation;

    my $contig_id = "";
    my $fset_id   = "";

    if (defined($contig->id)) {
       $contig_id = $contig->id;
    }
    if (defined($fset->id)) {
       $fset_id = $fset->id;
    }

    $translation->temporary_id($contig_id . "." . $fset_id);
    $translation->start        (1);
    $translation->end          ($exons[$#exons]->end - $exons[$#exons]->start + 1);
    $translation->start_Exon($exons[0]);
    $translation->end_Exon  ($exons[$#exons]);
    $transcript->translation($translation);

 #  print STDERR "Phase $startphase " . $transcript->translate->seq . "\n";

    $startphase++;
  }
  #print "finshed  fset2transcript_3frame\n";
  return @transcripts;
}


sub fset2transcript_with_seq {
    my ($genscan,$seq)=@_;

  
    unless ($genscan->isa ("Bio::EnsEMBL::SeqFeatureI"))
    {print "$genscan must be Bio::EnsEMBL::SeqFeatureI\n";}
    unless ($seq->isa ("Bio::PrimarySeqI") || $seq->isa ("Bio::SeqI"))
    {print "$seq must be Bio::SeqI or a Bio::PrimarySeqI\n";}

    #print STDERR "running fset2transcript\n";
    my $transcript = new Bio::EnsEMBL::Transcript;
    $transcript->temporary_id($seq->id . "." . $genscan->seqname);
        
    my @exons;
    my $count= 1;
    
    foreach my $f ($genscan->sub_SeqFeature) {
  
	my $exon  = new Bio::EnsEMBL::Exon;
        $exon->contig   ($seq);
	$exon->start    ($f->start);
	$exon->end      ($f->end  );
	$exon->strand   ($f->strand);
	$exon->phase    ($f->phase);
	$exon->end_phase( ($exon->phase + $exon->length)%3 );
	#$exon->score ($f->score);
	#print STDERR "contig is a = ".$seq."\n";
	$exon->slice($seq);
	
	push(@exons,$exon);
	$count++;
	
    }

    foreach my $exon (@exons) {
       	
      $transcript->add_Exon($exon);
	
	
    }
    return $transcript;
   
}



1;
