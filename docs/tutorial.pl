#!/usr/local/bin/perl -w

# This file is the companion script to the ensembl_tutorial and has been tested with
# the homo_sapiens_core_130 databases and the branch-ensembl-3 code.

# It tests the major parts of the ensembl API and is intended to be
# used as an introduction to new users.

# MC Jan 2002
# Dan Andrews Aug 2002

use Bio::EnsEMBL::Map::DBSQL::Obj;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::ExternalData::ESTSQL::DBAdaptor;
use Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor;
use Bio::EnsEMBL::TranscriptFactory;

use Bio::SeqIO;

use strict;

my $host   = 'kaka.sanger.ac.uk';
my $user   = 'anonymous';
my $host    = 'ecs1d',
my $user    = 'ensro';
my $dbname = 'current';
my $dbname = 'homo_sapiens_core_130';
my $mapname= 'homo_sapiens_maps_130';
my $snpname= 'homo_sapiens_snp_7_29';
my $estname= 'homo_sapiens_est_7_29';
my $path   = 'NCBI_29';

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $host,
					    -user   => $user,
					    -dbname => $dbname);


#my @clones = $db->get_all_Clone_id;                   ###########  Deprecated method with no obvious replacement


my $clone = $db->get_CloneAdaptor()->fetch_by_accession('AC005663');


print "Clone is " . $clone->id . "\n";


my @contigs = $clone->get_all_Contigs;


#foreach my $contig (@contigs) {
#
#  my $length = $contig->length;
#  my $id     = $contig->id;
#
#  print $contig->seq . "\n";
#  print $contig->subseq(1,100,1) . "\n";                    
#
#  my $seqio = new Bio::SeqIO(
#			      '-fh'     => \*STDOUT,
#			      '-format' => 'fasta');
#
#  $seqio->write_seq($contig);
#}

my $contig = $contigs[0];

#my $maskedseq = $contig->get_repeatmasked_seq;           ########## Doesn't appear masked - check against REAL data.

#print $maskedseq->seq . "\n";

#my @repeats = $contig->get_all_RepeatFeatures;
#
#foreach my $repeat (@repeats) {
#  print "Feature is " . $repeat->gffstring . "\n";
#}
#
#foreach my $repeat ($contig->get_all_RepeatFeatures) {
#
#   print "Name   : " . $repeat->seqname . "\n";
#   print "Start  : " . $repeat->start   . "\n";
#   print "End    : " . $repeat->end     . "\n";
#   print "Strand : " . $repeat->strand . "\n";
#   print "Score  : " . $repeat->score  . "\n";
#
#}


#foreach my $repeat ($contig->get_all_RepeatFeatures) {
#	  print "Hit start " . $repeat->hstart . "\n";
#	  print "Hit end   " . $repeat->hend . "\n";
#}

#foreach my $feat ($contig->get_all_SimilarityFeatures) {
#  # need to protect this for anything that is not a
#  # similaritypair, cpgs are in current db for similarity
#  next if( $feat->source_tag eq 'cpg' );
#  print "Feature scores are " . $feat->score      . "\t" .
#				  $feat->percent_id . "\t" .
#				  $feat->p_value    . "\n";
#
#  my $analysis = $feat->analysis;
#
#  print "Database :   " . $analysis->db . "\n";
#  print "Program :    " . $analysis->program . "\n";
#  print "Parameters : " . $analysis->parameters . "\n";
#
#}


#my $mapdb = new Bio::EnsEMBL::Map::DBSQL::Obj(-host   => $host,       ################# This is untested.
#					      -user   => $user,
#					      -dbname => $mapname);
#
#$db->mapdbname($mapname);
#$db->mapdb($mapdb);
#
#my @markers = $contig->get_landmark_MarkerFeatures;                ########### Presently not implemented in RawContig

#foreach my $marker (@markers) {
#  print $marker->gffstring . "\n";
#}
                                                                    ########### This bit up to the previous message is untested.
my $gene_adaptor = $db->get_GeneAdaptor;
my @genes = $gene_adaptor->fetch_by_contig_list(($contig->id));

#foreach my $gene (@genes) {
#  print "Gene : " . $gene->stable_id . "\n";
#}


#foreach my $gene (@genes) {
#  if ($gene->is_known) {
#    my @dblinks = $gene->each_DBLink;
#    
#    foreach my $link (@dblinks) {
#      print "Gene " . $gene->stable_id . " links to " . 
#	   $link->display_id . " " .
#	     $link->database . "\n";
#      
#      my @syns = $link->get_synonyms;
#     
#      print "Synonyms for gene are @syns\n" if scalar @syns > 0;
#    }
#  } else {
#    print "Gene " . $gene->stable_id . " is not a known gene\n";
#  }
#  my $description = $gene->description;
#
#  print "Gene description is $description\n" if defined $description;
#
#}
  
my $gene = $genes[0];

my @transcripts = $gene->get_all_Transcripts();
my $transcript = $transcripts[0];

#foreach my $exon ($transcript->get_all_Exons) {
#   print "Found exon " . $exon->stable_id    . " " .
#			   $exon->start . " " .
#			   $exon->end   . " " .
#			   $exon->seq->seq . "\n";
#}
#
#my $seqio = new Bio::SeqIO(
#			   '-format' => 'fasta',
#			   '-fh'     => \*STDOUT);

#foreach my $transcript ($gene->get_all_Transcripts) {
#  my $peptide = $transcript->translate;
#
#  $seqio->write_seq($peptide);
#}

#@genes = $contig->get_all_Genes('evidence');                  ########## This will be done in another way soon.
#
#foreach my $gene (@genes) {
#  foreach my $transcript ($gene->each_Transcript) {
#    foreach my $exon ($transcript->get_all_Exons) {
#      my @evidence = $exon->each_Supporting_Feature;
#
#      foreach my $f (@evidence) {
#	 print "Evidence " . $f->gffstring . "\n";
#      }
#    }
#  }
#}                                                                ##########  The above block of code is untested.

#my @predicted_transcripts = $contig->get_all_PredictionFeatures;
#
#foreach my $transcript (@predicted_transcripts) {
#    my @exons = $transcript->get_all_Exons;
#    print "Genscan prediction has " . scalar(@exons) . " exons\n";
#    foreach my $exon (@exons) {
#	   print $exon->start  . " - " .
#		 $exon->end    . " : " .
#		 $exon->strand . " " .
#		 $exon->phase ."\n";
#    }
#}
#
#if (scalar @predicted_transcripts > 0) {
#    print "Genscan Peptide is " . $predicted_transcripts[0]->translate . "\n";
#}

#$db->assembly_type($path);
#
#my $slice_adaptor = $db->get_SliceAdaptor;
#my $slice = $slice_adaptor->fetch_by_chr_start_end('1',1,100000);
#
#my $slice2 = $slice_adaptor->fetch_by_contig_name('AC005663.2.1.103122',10000);
#
#my $slice3 = $slice_adaptor->fetch_by_clone_accession('AC005663',1000);
#
#my $slice4 = $slice_adaptor->fetch_by_gene_stable_id('ENSG00000099889',5000);
#
#my @rept = $slice->get_all_RepeatFeatures('RepeatMask');
#
#my @pred = $slice->get_all_PredictionTranscripts;
#
#my @sims = $slice->get_all_SimilarityFeatures_above_score('Vertrna', 0.00001); 
#
#my @slice_genes = $slice->get_all_Genes;
#
#my $chrname  = $slice->chr_name;
#my $chrstart = $slice->chr_start;
#my $chrend   = $slice->chr_end;
#
#print "Chromosome " . $chrname . "  Start/end " . $chrstart . " " . $chrend . "\n";


my $translation = $transcript->translation;

#$translation->start_exon->stable_id;
#$translation->end_exon->stable_id;
#
#$translation->start;
#$translation->end;
#
#print "Translation " . $translation->stable_id . ":" .$translation->dbID . " " . $translation->start . " " . $translation->end . "\n";

#print "Fetching interpro\n";
#
#my @interpro = $db->get_GeneAdaptor->get_Interpro_by_geneid('ENSG00000099889');
#
#
#my $protein_adaptor = $db->get_ProteinAdaptor();
#
##my $protein  = $protein_adaptor->fetch_Protein_by_dbid($translation->dbID);
#
#my $protein = $protein_adaptor->fetch_by_Transcript_id($transcript->stable_id);
#
#
#my @prot_features = $protein->all_SeqFeature;
#
#foreach my $pf (@prot_features) {
#  print $pf->gffstring . "\n";
#}

my $host2 = 'kaka.sanger.ac.uk';
my $user2 = 'anonymous';

my $est= new Bio::EnsEMBL::ExternalData::ESTSQL::DBAdaptor
  (-host   => $host2,          ########## temporary setting
   -dbname => $estname,
   -user   => $user2);          ########## temporary setting

my $est_external_feature_factory = $est->get_EstAdaptor;          ########## DBI errors due to using old db with new schema
                                                                   ########## test with a new db built with the new schema.
$db->add_ExternalFeatureFactory($est_external_feature_factory);

my $snpdb = new Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor
  (-host   => $host2,          ########## temporary setting
   -dbname => $snpname,
   -user   => $user2);          ########## temporary setting

$db->add_ExternalFeatureFactory($snpdb);

my @ext_features = $contig->get_all_ExternalFeatures;

foreach my $ext (@ext_features) {
  if ($ext->isa("Bio::EnsEMBL::FeaturePair")) {
    print "FP " . $ext->gffstring . "\n";
  } elsif ($ext->isa("Bio::EnsEMBL::ExternalData::Variation")) {
    print "VAR " . $ext->id . "\t" . $ext->status . "\t".  $ext->clone_name . "\t" . $ext->start_in_clone_coord . "\n";
  }
}











