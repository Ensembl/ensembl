#!/usr/local/bin/perl

# This file is the companion script to the ensembl_tutorial and has been tested with
# the homo_sapiens_core_130 databases and the branch-ensembl-3 code.

# It tests the major parts of the ensembl API and is intended to be
# used as an introduction to new users.

# MC Jan 2002

use Bio::EnsEMBL::Map::DBSQL::Obj;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::ExternalData::ESTSQL::DBAdaptor;
use Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor;

use Bio::SeqIO;

use strict;

my $host   = 'kaka.sanger.ac.uk';
my $user   = 'anonymous';
#my $host    = 'ecs1d',
#my $user    = 'ensro';
my $dbname = 'homo_sapiens_core_7_29';
my $mapname= 'homo_sapiens_maps_130';
my $snpname= 'homo_sapiens_snp_130';
my $estname= 'homo_sapiens_est_130';
my $path   = 'NCBI_29';

my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $host,
					    -user   => $user,
					    -dbname => $dbname);


#my @clones = $db->get_all_Clone_id;


my $clone = $db->get_Clone('AC005663');


print "Clone is " . $clone->id . "\n";


my @contigs = $clone->get_all_Contigs;


foreach my $contig (@contigs) {

  my $seqobj = $contig->primary_seq;
  my $length = $contig->length;
  my $id     = $contig->id;

  print $seqobj->seq . "\n";
  print $seqobj->subseq(1,100) . "\n";

  my $seqio = new Bio::SeqIO(-fh     => \*STDOUT,
			     -format => 'fasta');

  $seqio->write_seq($seqobj);;
}

my $contig = $contigs[0];

my $maskedseq = $contig->get_repeatmasked_seq;

my $contig = $contigs[0];

my @repeats = $contig->get_all_RepeatFeatures;


foreach my $repeat (@repeats) {
  print "Feature is " . $repeat->gffstring . "\n";
}




foreach my $repeat ($contig->get_all_RepeatFeatures) {

   print "Name   : " . $repeat->seqname . "\n";
   print "Start  : " . $repeat->start   . "\n";
   print "End    : " . $repeat->end     . "\n";
   print "Strand : " . $repeat->strand . "\n";
   print "Score  : " . $repeat->score  . "\n";

}



foreach my $repeat ($contig->get_all_RepeatFeatures) {
        print "Hit name  " . $repeat->hseqname . "\n";
        print "Hit start " . $repeat->hstart . "\n";
        print "Hit end   " . $repeat->hend . "\n";
}

foreach my $feat ($contig->get_all_SimilarityFeatures) {
  # need to protect this for anything that is not a
  # similaritypair, cpgs are in current db for similarity
  next if( $feat->source_tag eq 'cpg' );
  print "Feature scores are " . $feat->score      . "\t" .
                                $feat->percent_id . "\t" .
                                $feat->p_value    . "\n";

  my $analysis = $feat->analysis;

  print "Database :   " . $analysis->db . "\n";
  print "Program :    " . $analysis->program . "\n";
  print "Parameters : " . $analysis->parameters . "\n";

}

my $mapdb = new Bio::EnsEMBL::Map::DBSQL::Obj(-host   => $host,
					      -user   => $user,
					      -dbname => $mapname);

$db->mapdbname($mapname);
$db->mapdb($mapdb);

my @markers = $contig->get_landmark_MarkerFeatures;

foreach my $marker (@markers) {
  print $marker->gffstring . "\n";
}

my @genes = $contig->get_all_Genes;

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
      
      print "Synonyms for gene are @syns\n";
    }
  } else {
    print "Gene " . $gene->stable_id . " is not a known gene\n";
  }
  my $description = $gene->description;

  print "Gene description is $description\n";

}
  
my $gene = $genes[0];



my @transcripts = $gene->each_Transcript;
my $transcript = @transcripts[0];

foreach my $exon ($transcript->get_all_Exons) {
   print "Found exon " . $exon->stable_id    . " " .
                         $exon->start . " " .
                         $exon->end   . " " .
                         $exon->seq->seq . "\n";
}

my $seqio = new Bio::SeqIO(-format => 'fasta',
                           -fh     => \*STDOUT);

foreach my $transcript ($gene->each_Transcript) {
  my $peptide = $transcript->translate;

  $seqio->write_seq($peptide);
}

my @genes = $contig->get_all_Genes('evidence');

foreach my $gene (@genes) {
  foreach my $transcript ($gene->each_Transcript) {
    foreach my $exon ($transcript->get_all_Exons) {
      my @evidence = $exon->each_Supporting_Feature;

      foreach my $f (@evidence) {
	print "Evidence " . $f->gffstring . "\n";
      }
    }
  }
}

my @genscan = $contig->get_all_PredictionFeatures;


foreach my $genscan (@genscan) {
  my @exons = $genscan->sub_SeqFeature;
  print "Genscan prediction has " . scalar(@exons) . " exons\n";
  foreach my $exon (@exons) {
    print $exon->start  . " - " .
          $exon->end    . " : " .
          $exon->strand . " " .
          $exon->phase ."\n";
  }
}
print "Making transcript\n";
my $gen_transcript = Bio::EnsEMBL::DBSQL::Utils::fset2transcript($genscan[1],$contig);

print "Genscan Peptide is " . $gen_transcript->translate->seq . "\n";


$db->static_golden_path_type($path);


my $sgp = $db->get_StaticGoldenPathAdaptor;


my $vcontig = $sgp->fetch_VirtualContig_of_contig('AC005663.2.1.103122',10000);
print "Got contig\n";

my $vcontig = $sgp->fetch_VirtualContig_of_clone('AC005663',1000);
print "Got clone\n";

my $vcontig = $sgp->fetch_VirtualContig_of_gene('ENSG00000099889',5000);

print "Got gene\n";

my @rept  = $vcontig->get_all_RepeatFeatures;

print "Repeats " . scalar(@rept) . "\n";

my @pred  = $vcontig->get_all_PredictionFeatures;

print "Predictions " . scalar(@pred) . "\n";

my @sims  = $vcontig->get_all_SimilarityFeatures;

print "Similarities " . scalar(@sims) . "\n";

my @genes = $vcontig->get_all_Genes_exononly;

print "Genes " . scalar(@genes) . "\n";

my @contigs  = $vcontig->_vmap->each_MapContig;

print "Contigs " . scalar(@contigs) . "\n";

my $chrname  = $vcontig->_chr_name;
my $chrstart = $vcontig->_global_start;
my $chrend   = $vcontig->_global_end;

print "Start/end " . $chrname . " " . $chrstart . " " . $chrend . "\n";


my $translation = $transcript->translation;

$translation->start_exon->stable_id;
$translation->end_exon->stable_id;

$translation->start;
$translation->end;

print "Translation " . $translation->stable_id . ":" .$translation->dbID . " " . $translation->start . " " . $translation->end . "\n";

print "Fetching interpro\n";

my @interpro = $db->get_GeneAdaptor->get_Interpro_by_geneid('ENSG00000099889');


my $protein_adaptor = $db->get_Protein_Adaptor($db);

my $protein  = $protein_adaptor->fetch_Protein_by_dbid($translation->dbID);

my $protein = $protein_adaptor->fetch_Protein_by_transcriptId($transcript->stable_id);


my @prot_features = $protein->all_SeqFeature;

foreach my $pf (@prot_features) {
  print $pf->gffstring . "\n";
}

my $est= new Bio::EnsEMBL::ExternalData::ESTSQL::DBAdaptor
  (-host   => $host,
   -dbname => $estname,
   -user   => $user);

my $est_external_feature_factory = $est->get_EstAdaptor;

$db->add_ExternalFeatureFactory($est_external_feature_factory);

my $snpdb = new Bio::EnsEMBL::ExternalData::SNPSQL::DBAdaptor
  (-host   => $host,
   -dbname => $snpname,
   -user   => $user);

$db->add_ExternalFeatureFactory($snpdb);

my @ext_features = $contig->get_all_ExternalFeatures;

foreach my $ext (@ext_features) {
  if ($ext->isa("Bio::EnsEMBL::FeaturePair")) {
    print $ext->gffstring . "\n";
  } elsif ($ext->isa("Bio::EnsEMBL::ExternalData::Variation")) {
    print $ext->id . "\t" . $ext->status . "\t".  $ext->clone_name . "\t" . $ext->start_in_clone_coord . "\n";
  }
}


