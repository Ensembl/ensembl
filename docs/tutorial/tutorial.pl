#!/usr/local/bin/perl -w

# This file is the companion script to the ensembl_tutorial and has been tested with
# the homo_sapiens_core_M_29 database on ecs1d and the main trunk code.

# It tests the major parts of the ensembl API and is intended to be
# used as an introduction to new users.

# MC Jan 2002
# Updated to main trunk code - Dan Andrews Oct 2002

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Lite::DBAdaptor;
use Bio::EnsEMBL::ExternalData::ESTSQL::DBAdaptor;
use Bio::SeqIO;

use strict;

### Currently there is no trunk db on kaka.
my $host = 'ecs1d';
my $user = 'ensro';
my $dbname = 'homo_sapiens_core_9_30';

### External dbs
my $estname  = 'homo_sapiens_est_M_29';
my $litename = 'homo_sapiens_lite_M_29';


my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $host,
					    -user   => $user,
					    -dbname => $dbname);

my $clone = $db->get_CloneAdaptor()->fetch_by_accession('AC005663');

print "Clone is " . $clone->id . "\n";

my @contigs = $clone->get_all_Contigs;

foreach my $contig (@contigs) {

  my $length = $contig->length;
  my $id     = $contig->id;

  print $contig->seq . "\n";
  print $contig->subseq(1,100,1) . "\n";                    

  my $seqio = new Bio::SeqIO(
			      '-fh'     => \*STDOUT,
			      '-format' => 'fasta');

  $seqio->write_seq($contig);
}

my $contig = $contigs[0];

my $maskedseq = $contig->get_repeatmasked_seq;
print $maskedseq->seq . "\n";

#$db->assembly_type($path);

my $slice_adaptor = $db->get_SliceAdaptor;
my $slice = $slice_adaptor->fetch_by_chr_start_end('1',1,100000);

my $slice2 = $slice_adaptor->fetch_by_contig_name('AC005663.2.1.103122',10000);

my $slice3 = $slice_adaptor->fetch_by_clone_accession('AC005663',1000);

my $slice4 = $slice_adaptor->fetch_by_gene_stable_id('ENSG00000099889',5000);

my @rept = $slice->get_all_RepeatFeatures('RepeatMask');

my @pred = $slice->get_all_PredictionTranscripts;

my $chrname  = $slice->chr_name;
my $chrstart = $slice->chr_start;
my $chrend   = $slice->chr_end;

print "Chromosome " . $chrname . "  Start/End " . $chrstart . " " . $chrend . "\n";

my @repeats = $slice->get_all_RepeatFeatures;

foreach my $repeat (@repeats) {
  print $repeat->gffstring . "\n";
}

foreach my $repeat ($slice->get_all_RepeatFeatures) {

   print "Name   : " . $repeat->seqname . "\n";
   print "Start  : " . $repeat->start   . "\n";
   print "End    : " . $repeat->end     . "\n";
   print "Strand : " . $repeat->strand . "\n";
   print "Score  : " . $repeat->score  . "\n";

}


foreach my $repeat ($slice->get_all_RepeatFeatures) {
	  print "Hit start " . $repeat->hstart . "\n";
	  print "Hit end   " . $repeat->hend . "\n";
}

foreach my $feat ($slice->get_all_DnaAlignFeatures_above_score('Vertrna', 0.00001)) {
  # need to protect this for anything that is not a
  # similaritypair, cpgs are in current db for similarity
  next if( $feat->source_tag eq 'cpg' );
  print "Feature scores are " . $feat->score      . "\t" .
				    $feat->percent_id . "\t" .
				    $feat->p_value    . "\n";

  my $analysis = $feat->analysis;

  print "Analysis:    " . $analysis->gff_source . "\n";
  print "Database :   " . $analysis->db . "\n";
  print "Program :    " . $analysis->program . "\n";
  print "Parameters : " . $analysis->parameters . "\n";

}

my @prot_feats = $slice->get_all_ProteinAlignFeatures_above_score('Swall', 0.00001);

### This is commented out as it takes yonks to run.  It illustrates the point if necessary though.

#my @sims = $slice->get_all_DnaAlignFeatures_above_score('Vertrna',0.00001);
#my @reps = $slice->get_all_RepeatFeatures;

#foreach my $sim (@sims){
#  foreach my $rep (@reps){
#    if ($sim->overlaps($rep)){
#      print "Similarity feature found with an overlapping Repeat feature\n";
#    }
#  }
#}

my @genes = $slice->get_all_Genes;

foreach my $gene (@genes) {
  print "Gene : " . $gene->stable_id . "\n";
}


foreach my $gene (@genes) {
  if ($gene->is_known) {
    my @dblinks = $gene->each_DBLink;
  
    foreach my $link (@dblinks) {
      print "Gene " . $gene->stable_id . " links to " . 
	   $link->display_id . " " .
	     $link->database . "\n";
    
      my @syns = $link->get_synonyms;
    
      print "Synonyms for gene are @syns\n" if scalar @syns > 0;
    }
  } else {
    print "Gene " . $gene->stable_id . " is not a known gene\n";
  }
  my $description = $gene->description;

  print "Gene description is $description\n" if defined $description;
}
  
my $gene = $genes[0];

my @transcripts = $gene->get_all_Transcripts();
my $transcript = $transcripts[0];

foreach my $exon ($transcript->get_all_Exons) {

   print "Found exon " . $exon->stable_id    . " " .
			   $exon->start . " " .
			   $exon->end   . " " .
			   $exon->seq->seq . "\n";

  my @evidence = $exon->each_Supporting_Feature;

  foreach my $item (@evidence) {
    print "Evidence " . $item->gffstring . "\n";
  }
}

my $seqio = new Bio::SeqIO(
			   '-format' => 'fasta',
			   '-fh'     => \*STDOUT);

foreach my $transcript ($gene->get_all_Transcripts) {
  my $peptide = $transcript->translate;

  $seqio->write_seq($peptide);
}

#my @predicted_transcripts = $contig->get_all_PredictionFeatures;
my @predicted_transcripts = $slice->get_all_PredictionTranscripts;


foreach my $transcript (@predicted_transcripts) {
    my @exons = $transcript->get_all_Exons;
    print "Genscan prediction has " . scalar(@exons) . " exons\n";
    foreach my $exon (@exons) {
	   print $exon->start  . " - " .
		 $exon->end    . " : " .
		 $exon->strand . " " .
		 $exon->phase ."\n";
    }
    print "Genscan Peptide is " . $transcript->translate . "\n";
}


my $translation = $transcript->translation;

$translation->start_exon->stable_id;
$translation->end_exon->stable_id;

$translation->start;
$translation->end;

print "Translation " . $translation->stable_id . ":" .$translation->dbID . " " . $translation->start . " " . $translation->end . "\n";

print "Fetching interpro\n";

my @interpro = $db->get_GeneAdaptor->get_Interpro_by_geneid('ENSG00000099889');

my $protein_adaptor = $db->get_ProteinAdaptor();

#my $protein  = $protein_adaptor->fetch_by_translation_id($translation->dbID);
my $protein = $protein_adaptor->fetch_by_transcript_id($transcript->stable_id);

my @prot_features = $protein->all_SeqFeature;

foreach my $pf (@prot_features) {
  print $pf->gffstring . "\n";
}


my $est_db= Bio::EnsEMBL::ExternalData::ESTSQL::DBAdaptor->new(-host   => $host,
							    -dbname => $estname,
							    -user   => $user);


$db->add_db_adaptor('est', $est_db);
$est_db->add_db_adaptor('core', $db);

my @est_features = $slice->get_all_DnaAlignFeatures_above_score('ex_e2g_feat', 0.0001);

foreach my $est (@est_features) {
    print " " . $est->gffstring . "\n";
}


my $lite_db = Bio::EnsEMBL::Lite::DBAdaptor->new(-host   => $host,
						 -user   => $user,
						 -dbname => $litename);

$db->add_db_adaptor('lite', $lite_db);
$lite_db->add_db_adaptor('core', $db);

my @snp_features = $slice->get_all_SNPs;

foreach my $snp (@snp_features) {
    print "snp " . $snp->position . "\n";
}

my @landmark_features = $slice->get_landmark_MarkerFeatures;

foreach my $marker (@landmark_features) {
    print $marker->gffstring . "\n";
}


