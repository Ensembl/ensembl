#!/usr/local/ensembl/bin/perl -w

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Transcript;
use Bio::EnsEMBL::Gene;
use Bio::EnsEMBL::Translation;
use Bio::EnsEMBL::Exon;

$| = 1;

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => "ecs1d",
					    -user => "ecs1dadmin",
					    -dbname => "elegans_newschema",
					    -pass  => "TyhRv",
					   );
#CHROMOSOME_I	curated CDS	11641	11689	.	+	0	Sequence "Y74C9A.2"
#CHROMOSOME_I	curated CDS	14951	15160	.	+	2	Sequence "Y74C9A.2"
#CHROMOSOME_I	curated CDS	16473	16585	.	+	2	Sequence "Y74C9A.2"
#CHROMOSOME_I	curated CDS	43733	43961	.	+	0	Sequence "Y74C9A.1"
#CHROMOSOME_I	curated CDS	44030	44234	.	+	2	Sequence "Y74C9A.1"
#CHROMOSOME_I	curated CDS	44281	44328	.	+	1	Sequence "Y74C9A.1"
#CHROMOSOME_I	curated CDS	44521	44677	.	+	1	Sequence "Y74C9A.1"
#CHROMOSOME_I	curated CDS	49921	50016	.	+	0	Sequence "Y48G1C.4"
#CHROMOSOME_I	curated CDS	50815	51030	.	+	0	Sequence "Y48G1C.4"
#CHROMOSOME_I	curated CDS	52283	52410	.	+	0	Sequence "Y48G1C.4"
#CHROMOSOME_I	curated CDS	52466	52572	.	+	1	Sequence "Y48G1C.4"


my %chromosome;
my $sql = "select name, length from chromosome";
my $sth = $db->prepare($sql);
$sth->execute;

my $slice_adaptor = $db->get_SliceAdaptor();

while(my ($name, $length) = $sth->fetchrow){
  my $slice = $slice_adaptor->fetch_by_chr_start_end($name, 1, $length);
  $chromosome{$name} = $slice;
} 



my $analysis_adaptor = $db->get_AnalysisAdaptor();
my $analysis = $analysis_adaptor->fetch_by_logic_name('wormbase');
my $forward_count = 0;
my $reverse_count = 0;
my $count = 0;


my %genes;
my @stored;

my $file = shift;
open(FH, $file) or die "couldn't open $file $!";

my $transcripts = &process_file(\*FH);
my $processed_transcripts = &process_transcripts($transcripts, $analysis);

#print "there are ".keys(%$processed_transcripts)." transcript\n";

my $genes = &create_transcripts($processed_transcripts);

#print "there are ".keys(%$genes)." genes\n";


my $gene_adaptor = $db->get_GeneAdaptor;
my @genes;
foreach my $gene_id(keys(%$genes)){
  my $db_id;
  my $transcripts = $genes->{$gene_id};
  my $unpruned = &create_gene($transcripts, $gene_id);
##  print STDERR "gene ".$unpruned."\n";
  my $gene = &prune_Exons($unpruned);
  push(@genes, $gene);
  #eval{
  #  $gene->transform;
  #};
  foreach my $transcript (@{$gene->get_all_Transcripts}){
    print "transcript ".$transcript->stable_id."\n";
    my @exons = @{$transcript->get_all_Exons};
    if($exons[0]->strand == -1){
      @exons = sort{$b->start <=> $a->start} @exons;
    }else{
      @exons = sort{$a->start <=> $b->start} @exons;
    }
    foreach my $exon(@exons){
      print "exon ".$exon->stable_id."\t".$exon->seqname."\t ".$exon->start."\t ".$exon->end."\t ".$exon->strand."\t ".$exon->phase."\t ".$exon->end_phase."\n";
    }
    
  }

#  if($@){
#    warn "couldn't transform ".$gene." $@\n";
#    next;
#  }
#  eval{
#    #$gene_adaptor->store($gene);
#  };
#  if($@){
#    warn "couldn't store ".$gene."\n";
#    next;
#  }
#  $db_id = $gene->dbID;
#  print STDERR "have stored gene ".$gene->dbID."\n";
#  #my $stored_gene = $gene_adaptor->fetch_by_dbID($db_id);
#  #&translation_check($stored_gene);
}
 
print "there are ".@genes." genes in wormbase 90\n";

sub process_file{
  my ($fh) = @_;
  
  my %transcripts;
 LOOP: while(<$fh>){
#    print;
    chomp;
    my($chr, $status, $type, $start, $end, $score, $strand, $frame, $sequence, $gene) = split;
    my $element = $_;
    if(!$status && !$type){
      #print "status and type no defined skipping\n";
      next LOOP;
    }
    my $line = $status." ".$type;
#    print "line ".$line."\n";
    if($line ne 'curated CDS'){
      next LOOP;
    }
    $gene =~ s/\"//g;
    if(!$transcripts{$gene}){
      $transcripts{$gene} = [];
      push(@{$transcripts{$gene}}, $element);
    }else{
      push(@{$transcripts{$gene}}, $element);
    }
    
  }
  return \%transcripts;
}


sub process_transcripts{
  my ($transcripts, $analysis) = @_;
  
  my %genes;
  my %transcripts = %$transcripts;
  my @names = keys(%transcripts);
  #print STDERR "PROCESSING TRANSCRIPTS \n";
  foreach my $name(@names){
    my @lines = @{$transcripts{$name}};
    $transcripts{$name} = [];
    my @exons;
    foreach my $line(@lines){
      my($chr, $status, $type, $start, $end, $score, $strand, $frame, $sequence, $gene) = split /\s+/, $line;
      $chr =~ s/CHROMOSOME_//;
      if($start == $end){
	next;
      }
      my $slice = $chromosome{$chr};
      my $exon = new Bio::EnsEMBL::Exon;
      my $phase = (3 - $frame)%3;
      #print "have ".$phase." phase\n";
      $exon->start($start);
      $exon->end($end);
      $exon->analysis($analysis);
      $exon->contig($slice);
      $exon->phase($phase);
      my $end_phase = ($phase + ($exon->end-$exon->start) + 1)%3;
      $exon->end_phase($end_phase);
      if($strand eq '+'){
	$exon->strand(1);
      }else{
	$exon->strand(-1);
      }
      $exon->score(100);
      push(@exons, $exon);
    }
    if($exons[0]->strand == -1){
      @exons = sort{$b->start <=> $a->start} @exons;
    }else{
      @exons = sort{$a->start <=> $b->start} @exons;
    }
    #print STDERR $name." transcript\n";
    my $phase = 0;
    foreach my $e(@exons){
      #print "phase ".$phase."\n";
      #$e->phase($phase);
      #my $end_phase = ($phase + ($e->end-$e->start) + 1)%3;
      #$e->end_phase($end_phase);
      #print "end phase ".$end_phase."\n";
      #$phase = $end_phase;
      push(@{$transcripts{$name}}, $e);
    }
  }
  #print STDERR "FINISHED PROCESSING TRANSCRIPTS\n";
  return \%transcripts;

}





sub create_transcripts{
  my ($transcripts) = @_;
  #print STDERR "CREATING TRANSCRIPTS\n";
  my %transcripts = %$transcripts;
  my @non_translate;
  my %genes;
  my $gene_name;
  my $transcript_id;
  foreach my $transcript(keys(%transcripts)){
    my $time = time;
    #print STDERR "transcript ".$transcript." \n";
    my @exons = @{$transcripts{$transcript}};
    if($transcript =~ /\w+\.\d+\w+/){
     ($gene_name) = $transcript =~ /(\w+\.\d+)\w+/;
     $transcript_id = $transcript;
    }else{
      $gene_name = $transcript;
      $transcript_id = $transcript;
    }
  
    my $transcript = new Bio::EnsEMBL::Transcript;
    my $translation = new Bio::EnsEMBL::Translation;
    my @sorted_exons;
    if($exons[0]->strand == 1){
      @sorted_exons = sort{$a->start <=> $b->start} @exons
    }else{
      @sorted_exons = sort{$b->start <=> $a->start} @exons  
    }
    my $exon_count = 1;
    my $phase = 0;
    foreach my $exon(@sorted_exons){
      $exon->created($time);
      $exon->modified($time);
      $exon->version(1);
      $exon->stable_id($transcript_id.".".$exon_count);
      #$exon->phase($phase);
      #my $end_phase = ($phase + ($end-$start) + 1)%3; 
      #$exon->end_phase($end_phase);
      #$phase = $end_phase;
      $exon_count++;
      $transcript->add_Exon($exon);
    }
    $translation->start_Exon($sorted_exons[0]);
    $translation->end_Exon  ($sorted_exons[$#sorted_exons]);
 
    if ($sorted_exons[0]->phase == 0) {
      $translation->start(1);
    } elsif ($sorted_exons[0]->phase == 1) {
      $translation->start(3);
    } elsif ($sorted_exons[0]->phase == 2) {
      $translation->start(2);
    }
    $translation->end  ($sorted_exons[$#sorted_exons]->end - $sorted_exons[$#sorted_exons]->start + 1);
#    $translation->version(1);
    $translation->stable_id($transcript_id);
    $transcript->translation($translation);
    $transcript->version(1);
    $transcript->stable_id($transcript_id);
   # my $mrna = $transcript->translateable_seq;
#    my $last_codon = substr($mrna, -3);
#    if($last_codon eq 'TAG' || $last_codon eq 'TGA' || $last_codon eq 'TAA'){
#      my $translation_end =   ($sorted_exons[$#sorted_exons]->end - $sorted_exons[$#sorted_exons]->start + 1)-3;
#      if($translation_end <= 0){
#	my $exon_number = $#sorted_exons
#      }
#      $translation->end  (($sorted_exons[$#sorted_exons]->end - $sorted_exons[$#sorted_exons]->start + 1)-3); 
#    }
    my $peptide = $transcript->translate->seq;
    if($peptide =~ /\*./){
     # warn "transcript ".$transcript->stable_id." doesn't translate and won't be used\n";
     # push(@non_translate, $transcript);
     # next;
    }
    #print STDERR "transcript ".$transcript->stable_id." has ".@{$transcript->get_all_Exons}." exons and produces and translation: ".$peptide."\n";
    my @stored_exons = @{$transcript->get_all_Exons};
    #&display_exons(@stored_exons);
    if(!$genes{$gene_name}){
      $genes{$gene_name} = [];
      push(@{$genes{$gene_name}}, $transcript);
    }else{
      push(@{$genes{$gene_name}}, $transcript);
    }
  }

  foreach my $t(@non_translate){
    print STDERR "transcript ".$t->stable_id." has ".@{$t->get_all_Exons}." exons\n";
    my @exons = @{$t->get_all_Exons};
    &display_exons(@exons);
    &non_translate($t);
  }
  #print STDERR "FINISHED CREATEING TRANSCRIPTS\n";
  return \%genes;

}


sub create_gene{
  my ($transcripts, $name) = @_;
  #print STDERR "CREATING GENE\n";
  my $time = time;
  my $gene = new Bio::EnsEMBL::Gene; 
  my $exons = $transcripts->[0]->get_all_Exons;
  my $analysis = $exons->[0]->analysis;
  $gene->analysis($analysis);
  $gene->type($analysis->logic_name);
  $gene->created($time);
  $gene->modified($time);
  $gene->version(1);
  $gene->stable_id($name);
  foreach my $transcript(@$transcripts){
    $gene->add_Transcript($transcript);
  }
  
  return $gene;
}


sub prune_Exons {
  my ($gene) = @_;
  
  my @unique_Exons; 
  
  # keep track of all unique exons found so far to avoid making duplicates
  # need to be very careful about translation->start_Exon and translation->end_Exon
  
  foreach my $tran (@{$gene->get_all_Transcripts}) {
    my @newexons;
    foreach my $exon (@{$tran->get_all_Exons}) {
      my $found;
      #always empty
    UNI:foreach my $uni (@unique_Exons) {
	if ($uni->start  == $exon->start  &&
	    $uni->end    == $exon->end    &&
	    $uni->strand == $exon->strand &&
	    $uni->phase  == $exon->phase  &&
	    $uni->end_phase == $exon->end_phase
	   ) {
	  $found = $uni;
	  last UNI;
	}
      }
      if (defined($found)) {
	push(@newexons,$found);
	if ($exon == $tran->translation->start_Exon){
	  $tran->translation->start_Exon($found);
	}
	if ($exon == $tran->translation->end_Exon){
	  $tran->translation->end_Exon($found);
	}
      } else {
	push(@newexons,$exon);
	push(@unique_Exons, $exon);
      }
    }          
    $tran->flush_Exon;
    foreach my $exon (@newexons) {
      $tran->add_Exon($exon);
    }
  }
  return $gene;
}

sub translation_check{
  my ($gene) = @_;

  my @transcripts = @{$gene->get_all_Transcripts};
  foreach my $t(@transcripts){
    my $pep = $t->translate->seq;
    if($pep =~ /\*/){
      print STDERR "transcript ".$t->stable_id." doesn't translate\n";
    }
  }
  
}


sub display_exons{
  my (@exons) = @_;

  @exons = sort{$a->start <=> $b->start || $a->end <=> $b->end} @exons;

  if($exons[0]->can ('hseqname')){
    foreach my $e(@exons){
      print $e->hseqname."\t ".$e->start."\t ".$e->end."\t ".$e->strand."\t ".$e->phase."\t ".$e->end_phase."\n";
    }
  }else{
     foreach my $e(@exons){
      print $e->seqname."\t ".$e->start."\t ".$e->end."\t ".$e->strand."\t ".$e->phase."\t ".$e->end_phase."\n";
    }
  }
}

sub non_translate{
  my (@transcripts) = @_;
  
  foreach my $t(@transcripts){
    
    my @exons = @{$t->get_all_Exons};
#    print "transcript sequence :\n".$t->seq."\n";
    foreach my $e(@exons){
      my $seq = $e->seq;
      my $pep0 = $seq->translate('*', 'X', 0);
      my $pep1 = $seq->translate('*', 'X', 1);
      my $pep2 = $seq->translate('*', 'X', 2);
      print "exon sequence :\n".$e->seq->seq."\n\n";
      print $e->seqname." ".$e->start." : ".$e->end." translation in 0 frame\n ".$pep0->seq."\n\n";
      print $e->seqname." ".$e->start." : ".$e->end." translation in 1 phase\n ".$pep2->seq."\n\n";
      print $e->seqname." ".$e->start." : ".$e->end." translation in 2 phase\n ".$pep1->seq."\n\n";
      print "\n\n";
      
    }
    
  }
}
