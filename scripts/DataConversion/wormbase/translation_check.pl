#!/usr/local/ensembl/bin/perl -w

use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;


$| = 1;

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host => "ecs1d",
					    -user => "ecs1dadmin",
					    -dbname => "elegans_newschema",
					    -pass  => "TyhRv",
					   );

my %chromosome;
my $sql = "select name, length from chromosome";
my $sth = $db->prepare($sql);
$sth->execute;

my $slice_adaptor = $db->get_SliceAdaptor();

while(my ($name, $length) = $sth->fetchrow){
  my $slice = $slice_adaptor->fetch_by_chr_start_end($name, 1, $length);
  $chromosome{$name} = $slice;
} 


foreach my $name(keys(%chromosome)){

my @genes;

my $slice = $chromosome{$name};
#$sql = 'select gene_id from gene';

#$sth = $db->prepare($sql);

#$sth->execute;

#while(my ($gene_id) = $sth->fetchrow){
#  my $gene = $db->get_GeneAdaptor->fetch_by_dbID($gene_id);
#  push(@genes, $gene);
  
#}

@genes = @{$slice->get_all_Genes};
#print STDERR "there are ".@genes." genes\n";



my $sticky_count = 0;
my $non_translate = 0;
foreach my $gene(@genes){

  #$sticky_count += &sticky_check($gene);
  $non_translate += &translation_check($gene);
  #my ($transcript) = @{$gene->get_all_Transcripts};
  #if($translate){
  #  my @exons = @{$transcript->get_all_Exons};
    #&display_exons(@exons);
    #&non_translate($transcript);
  #}
}

#print "there are ".$sticky_count." transcripts containing sticky_exons\n";
print "there are ".$non_translate." non translating transcripts in ".$name."\n";

}



sub sticky_check{
  my ($gene) = @_;
  
  my $count = 0;
  my $sticky = 0;
  my @transcripts = @{$gene->get_all_Transcripts};
  foreach my $t(@transcripts){
    foreach my $e(@{$t->get_all_Exons}){
      if($e->isa("Bio::EnsEMBL::StickyExon")){
	$sticky = 1;
      }
    }
    if($sticky){
      $count++;
    }
    $sticky = 0;
  }
  return $count;
}
sub translation_check{
  my ($gene) = @_;

  my $count = 0;
  my @transcripts = @{$gene->get_all_Transcripts};
  foreach my $t(@transcripts){
    my $pep = $t->translate->seq;
    if($pep =~ /\*/){
      print STDERR "\ntranscript ".$t->stable_id." doesn't translate: \n\n";
      $count++;
    }
  }
  return $count;
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
