## Bioperl Test Harness Script for Modules
##
# CVS Version
# $Id$


# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.t'

#-----------------------------------------------------------------------
## perl test harness expects the following output syntax only!
## 1..3
## ok 1  [not ok 1 (if test fails)]
## 2..3
## ok 2  [not ok 2 (if test fails)]
## 3..3
## ok 3  [not ok 3 (if test fails)]
##
## etc. etc. etc. (continue on for each tested function in the .t file)
#-----------------------------------------------------------------------

#
#NB: Are not tested the following methods: get_Protein_annseq, write_all_Protein_features, write_Protein_feature.
#

## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..13\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}


use lib 't';
use EnsTestDB;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::RawContig;
use Bio::EnsEMBL::DBSQL::FeatureAdaptor;
use Bio::EnsEMBL::SeqFeature;
use Data::Dumper;

my $counter = 1;
$loaded = 1;
my $debug = 0;
print "ok ", $counter++, "\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new();
    
# Load some data into the db
$ens_test->do_sql_file("t/featureAdaptor.dump");
    
# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;
print "ok " , $counter++ , "\n";

#Get a feature adaptor
my $fa=Bio::EnsEMBL::DBSQL::FeatureAdaptor->new($db);
print "ok ", $counter++, "\n";


#Test fetch_PredictionFeature_by_id
my $pf = $fa->fetch_PredictionFeature_by_id("194643");
$debug && print "fetch_PredictionFeature_by_id\t";
($pf->isa("Bio::EnsEMBL::SeqFeatureI")) ? (print "ok " . $counter++ . "\n") : (print "not ok " . $counter++ . "\n");
my $contig = $db->get_Contig($pf->seqname); 


#Test fetch_RepeatFeatures_by_RawContig
my @rf1 = $fa->fetch_RepeatFeatures_by_RawContig($contig->id);
$debug && print "fetch_RepeatFeatures_by_RawContig\t";
($rf1[0]->isa("Bio::EnsEMBL::Repeat")) ? (print "ok " . $counter++ . "\n") : (print "not ok " . $counter++ . "\n");
#print Data::Dumper->Dump([@rf1], ['*rf1']);


# Test fetch_by_hid
my @f1 = $fa->fetch_by_hid('IL5_HUMAN');
#print Data::Dumper->Dump([@f1], ['*f1']);
$debug && print "fetch_by_hid\t";
($f1[0]->isa("Bio::EnsEMBL::SeqFeatureI")) ? (print "ok " . $counter++ . "\n") : (print "not ok " . $counter++ . "\n");


# Test fetch_PredictionFeature_as_Transcript
my $tr = $fa->fetch_PredictionFeature_as_Transcript("194643");
#print Data::Dumper->Dump([$tr], ['tr']);
$debug && print "fetch_PredictionFeature_as_Transcript\t";
($tr->isa("Bio::EnsEMBL::Transcript")) ? (print "ok " . $counter++ . "\n") : (print "not ok " . $counter++ . "\n");


# Test delete_by_RawContig
$fa->delete_by_RawContig($contig);
eval { $tr = $fa->fetch_PredictionFeature_as_Transcript("194643"); };
$debug && print "delete_by_RawContig\t";
($@) ? (print "ok " . $counter++ . "\n") : (print "not ok " . $counter++ . "\n");


# Test _store_PredictionFeature
$fa->store($contig, $pf);
$debug && print "_store_PredictionFeature\t";
my @f2 = eval{$fa->fetch_by_hid("Internal Exon");};
($f2[0]->isa("Bio::EnsEMBL::SeqFeatureI")) ? (print "ok " . $counter++ . "\n") : (print "not ok " . $counter++ . "\n");


# Test _store_FeaturePair
$fa->store($contig, @f1);
$debug && print "_store_FeaturePair\t";
my @f3 = $fa->fetch_by_hid('IL5_HUMAN');
(scalar(@f1) == scalar(@f3)) ? (print "ok " . $counter++ . "\n") : (print "not ok " . $counter++ . "\n");


# Test delete_by_RawContig_id
$fa->delete_by_RawContig_id($contig->id);
$debug && print "delete_by_RawContig_id\t";
my @f4 = eval { $fa->fetch_by_hid('IL5_HUMAN'); };
(!@f4) ? (print "ok " . $counter++ . "\n") : (print "not ok " . $counter++ . "\n");


# Test _store_single_feature
my $sf = Bio::EnsEMBL::SeqFeature->new(-seqname => 'moos', -start => 1, -end => 100, -strand => 1, -score => 1000, -primary_tag => 'test', -source_tag => 'test', -analysis => $pf->analysis);
$fa->store($contig, $sf);
$debug && print "_store_single_feature\t";
@f3 = $fa->fetch_by_hid('__NONE__');
($f3[0]->isa("Bio::EnsEMBL::SeqFeatureI")) ? (print "ok " . $counter++ . "\n") : (print "not ok " . $counter++ . "\n");


#Test _store_Repeat
$fa->store($contig, @rf1);
my @rf2 = $fa->fetch_RepeatFeatures_by_RawContig($contig->id);
$debug && print "fetch_RepeatFeatures_by_RawContig\t";
(scalar(@rf1) == scalar(@rf2)) ? (print "ok " . $counter++ . "\n") : (print "not ok " . $counter++ . "\n");

