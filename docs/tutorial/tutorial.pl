#!/usr/local/bin/perl -w

# This file is the companion script to the ensembl_tutorial and has been 
# tested with the homo_sapiens_core_9_30 database on ecs1d and the main trunk 
# code.

# It tests the major parts of the ensembl API and is intended to be
# used as an introduction to new users.

# MC Jan 2002
# Updated to main trunk code - Dan Andrews Oct 2002
# Re-updated to main trunk code - Graham McVicker Oct 2002

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Lite::DBAdaptor;
use Bio::EnsEMBL::ExternalData::ESTSQL::DBAdaptor;
use Bio::SeqIO;

use strict;

### Currently there is no trunk db on kaka.
my $host = 'kaka.sanger.ac.uk';
my $user = 'anonymous';
my $dbname = 'homo_sapiens_core_9_30';

### External dbs
my $estname  = 'homo_sapiens_est_9_30';
my $litename = 'homo_sapiens_lite_9_30';


#connect to the core database
my $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(-host   => $host,
					    -user   => $user,
					    -dbname => $dbname);

#get some ObjectAdaptors from the core database
my $gene_adaptor = $db->get_GeneAdaptor;
my $clone_adaptor = $db->get_CloneAdaptor;

#retrieve a clone via its EMBL accession
my $clone = $clone_adaptor->fetch_by_accession('AC005663');

print "Clone is " . $clone->id . "\n";

#get the contigs from a clone
my $contigs = $clone->get_all_Contigs;

foreach my $contig (@$contigs) {
  my $length = $contig->length;
  my $id     = $contig->id;

  #the sequence of a contig is obtained easily:
  print $contig->seq . "\n";

  #we can get a substring of this sequence too:
  print $contig->subseq(1,100,1) . "\n";                    

  #if we want to  write the sequence to a file we use bioperl again
  # (Note: if we are creating a new object we need to include a use to
  # Bio::SeqIO line at the start of our file)

  my $seqio = new Bio::SeqIO('-fh'     => \*STDOUT,
			     '-format' => 'fasta');

  $seqio->write_seq($contig);
}

my $contig = $contigs->[0];

#get a repeat masked sequence from a contig
my $maskedseq = $contig->get_repeatmasked_seq;
print $maskedseq->seq . "\n";

#we can create slices in several different ways:
my $slice_adaptor = $db->get_SliceAdaptor;

#create a slice on a region of the assembly
my $slice = $slice_adaptor->fetch_by_chr_start_end('1',1,100000);

#create a slice surrounding a contig
my $slice2 = $slice_adaptor->fetch_by_contig_name('AC005663.2.1.103122',10000);

#create a slice surrounding a clone
my $slice3 = $slice_adaptor->fetch_by_clone_accession('AC005663',1000);

#create a slice around a gene
my $slice4 = $slice_adaptor->fetch_by_gene_stable_id('ENSG00000099889',5000);

#get the sequence string from a slice
my $sequence = $slice->seq;

#get some information about where the slice lies on the assembly:
my $chr_name  = $slice->chr_name;
my $chr_start = $slice->chr_start;
my $chr_end   = $slice->chr_end;

print "Chromosome $chr_name Start/End $chr_start $chr_end \n";

#get information about the tiling path underneath the slice
my $tiling_path = $slice->get_tiling_path;

foreach my $tile (@$tiling_path) {
  my $tile_slice = $tile->assembled_Seq->name;
  my $tile_ctg = $tile->component_Seq->name;

  my $slice_start = $tile->assembled_start;
  my $slice_end   = $tile->assembled_end;
  my $ctg_start   = $tile->component_start;
  my $ctg_end     = $tile->component_end;
  my $ctg_ori     = $tile->component_ori;

  print "[$tile_slice] $slice_start - $slice_end => " .
        "[$tile_ctg] $ctg_start - $ctg_end ($ctg_ori)\n";
}

#get the repeat features off of a slice
my $repeats = $slice->get_all_RepeatFeatures;

#print info about the repeats in the form of a gff string
foreach my $repeat (@$repeats) {
  print $repeat->gffstring . "\n";
}

#print info about features in a different way:
foreach my $repeat (@$repeats) {
  print "Name  : " . $repeat->seqname . "\n";
  print "Start : " . $repeat->start   . "\n";
  print "End   : " . $repeat->end     . "\n";
  print "Strand: " . $repeat->strand  . "\n";
  print "Score : " . $repeat->score   . "\n";
}

#repeats also have coordinates on the hit sequence:
foreach my $repeat (@$repeats) {
  print "Hit start " . $repeat->hstart . "\n";
  print "Hit end   " . $repeat->hend   . "\n";
}

#get all of the similarity features off this slice:
my $features = $contig->get_all_SimilarityFeatures;

#get only wublastn blast hits against EMBL vertrna database
$features = $slice->get_all_DnaAlignFeatures('Vertrna', 0.00001);

foreach my $feat (@$features) {
  print "Feature scores are " . $feat->score      . "\t" .
                                $feat->percent_id . "\t" .
				$feat->p_value    . "\n";
}

#get only wublastp hits against SWALL database
$features = $slice->get_all_ProteinAlignFeatures('Swall', 0.0001);

#
#A bunch of different ways of getting features out of EnsEMBL follow:
#
my $simple_feature_adaptor = $db->get_SimpleFeatureAdaptor;

#get all of the eponine hits
my @eponine_hits = 
  @{$simple_feature_adaptor->fetch_all_by_RawContig($contig, 'Eponine')};

#get all cpg islands on a slice with score greater than 500
my $cpg_islands = 
  $simple_feature_adaptor->fetch_all_by_Slice_and_score($slice, 500, 'CpG');

#get the simple feature with the database identifier 1
my $simple_feature = $simple_feature_adaptor->fetch_by_dbID(1);

#get all of the simple features on a slice
my @simple_features = @{$slice->get_all_SimpleFeatures};

#get all of the tRNAscan hits on a contig
my $trnas = $contig->get_all_SimpleFeatures('tRNAscan');



#get the analysis information for a feature
my $feature = $features->[0];
my $analysis = $feature->analysis;
print "Analysis:    " . $analysis->gff_source . "\n";
print "Database :   " . $analysis->db . "\n";
print "Program :    " . $analysis->program . "\n";
print "Parameters : " . $analysis->parameters . "\n";


#get all the genes off of a slice:
my $genes = $slice->get_all_Genes;

foreach my $gene (@$genes) {
  print "Gene : " . $gene->stable_id . "\n";
}

# Get external descriptive information if a gene is known
foreach my $gene (@$genes) {
  if ($gene->is_known) {
    foreach my $link (@{$gene->get_all_DBLinks}) {
      print "Gene " . $gene->stable_id . " links to " . 
	  $link->display_id . " " .
	  $link->database . "\n";
      my @syns = @{$link->get_all_synonyms};
      print "Synonyms for gene are @syns\n" if scalar @syns > 0;
    }
  } else {
    print "Gene " . $gene->stable_id . " is not a known gene\n";
  }
}
  
my $gene = $genes->[0];

#genes also have a description
my $description = $gene->description;
print "Gene description is $description\n";

#get transcripts of a gene
my $transcripts = $gene->get_all_Transcripts();
my $transcript = $transcripts->[0];

#get exons of a transcript
foreach my $exon (@{$transcript->get_all_Exons}) {

   print "Found exon " . $exon->stable_id    . " " .
			   $exon->start . " " .
			   $exon->end   . " " .
			   $exon->seq->seq . "\n";
}


#dump the protein sequence of a transcript
my $seqio = new Bio::SeqIO('-format' => 'fasta',
			   '-fh'     => \*STDOUT);

foreach my $transcript (@{$gene->get_all_Transcripts}) {
  my $peptide = $transcript->translate;

  $seqio->write_seq($peptide);
}


#get the supporting evidence that was used to predict each exon
foreach my $gene (@$genes) {
  foreach my $trans (@{$gene->get_all_Transcripts}) {
    foreach my $exon (@{$trans->get_all_Exons}) {
      my $supporting_features = $exon->get_all_supporting_features;
      foreach my $feature (@$supporting_features) {
	print "Evidence " . $feature->gffstring . "\n";
      }
    }
  }
}


#get genscan predictions in the region of a slice
my $predicted_transcripts = $slice->get_all_PredictionTranscripts;

foreach my $ptrans (@$predicted_transcripts) {
    my $exons = $ptrans->get_all_Exons;
    print "Genscan prediction has " . scalar(@$exons) . " exon(s)\n";
    foreach my $exon (@$exons) {
	   print $exon->start  . " - " .
		 $exon->end    . " : " .
		 $exon->strand . " " .
		 $exon->phase ."\n";
    }
}

#print out a peptide predicted by genscan
my $ptrans = $predicted_transcripts->[0];
print "Genscan Peptide is " . $transcript->translate->seq . "\n";


#get the translation of a transcript
my $translation = $transcript->translation;

#the exon where translation starts
$translation->start_Exon->stable_id;
#the exon where translation ends
$translation->end_Exon->stable_id;

#the location within the start exon where translation starts
$translation->start;

#the location within the end exon where translation ends
$translation->end;

print "Translation " . $translation->stable_id . ":" .$translation->dbID . " " . $translation->start . " " . $translation->end . "\n";


my $protein_adaptor = $db->get_ProteinAdaptor();

#retrieve a protein via a translation id
my $protein  = $protein_adaptor->fetch_by_translation_id($translation->dbID);

#retrieve a protein via a transcript id
$protein = $protein_adaptor->fetch_by_transcript_id($transcript->stable_id);

#get all of the features of a protein
my $protein_features = $protein->get_all_ProteinFeatures;

foreach my $pf (@$protein_features) {
  print $pf->gffstring . "\n";
}

#connect to the EST database
my $est_db = 
  Bio::EnsEMBL::ExternalData::ESTSQL::DBAdaptor->new(-host   => $host,
						     -dbname => $estname,
						     -user   => $user);


$db->add_db_adaptor('est', $est_db);
$est_db->add_db_adaptor('core', $db);

#get EST alignment features off a slice
my $est_features = 
  $slice->get_all_DnaAlignFeatures('ex_e2g_feat');

foreach my $est (@$est_features) {
    print " " . $est->gffstring . "\n";
}


#connect to the lite database
my $lite_db = Bio::EnsEMBL::Lite::DBAdaptor->new(-host   => $host,
						 -user   => $user,
						 -dbname => $litename);

$db->add_db_adaptor('lite', $lite_db);
$lite_db->add_db_adaptor('core', $db);


# #print out snp features on a  slice
 my $snp_features = $slice->get_all_SNPs;

 foreach my $snp (@$snp_features) {
     print "snp " . $snp->start . "\n";
 }

#print out landmark marker features on a slice
my @landmark_features = @{$slice->get_all_landmark_MarkerFeatures};

foreach my $marker (@landmark_features) {
    print $marker->gffstring . "\n";
}


