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
BEGIN { $| = 1; print "1..19\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}


use lib 't';
use EnsTestDB;
use Bio::EnsEMBL::DBSQL::Feature_Obj;

$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new();
    
# Load some data into the db
$ens_test->do_sql_file("t/featobj.dump");
    
# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;
print "ok 2\n";

#Get a new feature_obj object
my $feature_obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($db);

#Get a contig object
#Before test method exeption
my $contig;
eval {
    my $clone=$db->get_Contig("wrong_id");
};
if ($@) {
    print "ok 3\n";
    
}
else {
    print "not ok 3\n";
    print STDERR "Trying to get a non existing clone not throw an exeption\n";
}
eval {
    $contig = $db->get_Contig("id_cont_test3");
};
if ($@) {
    print "not ok 4\n";
}
else {
    print "ok 4\n";
}

#Get an analysis object 
#Test the different values returned by the by the method
my $analysis;
eval {
    $analysis =$feature_obj->get_Analysis("wrong_id");
};
if ($@) {
    print "ok 5\n";
    
}
else {
    print "not ok 5\n";
    print STDERR "Trying to get a non existing analysis not throw an exeption\n";
}


eval {
    $analysis = $feature_obj->get_Analysis(1);
};
if ($@) {
    print "not ok 6\n";
    
}
else {
    print "ok 6\n";
}
if ($analysis->db eq "analysis_test3") {
    print "ok 7\n";
}
else {
    print "not ok7\n";
}

if ($analysis->db_version eq "nb3") {
    print "ok 8\n";
}
else {
    print "not ok8\n";
}

if ($analysis->program eq "gogene") {
    print "ok 9\n";
}
else {
    print "not ok 9\n";
}

if ($analysis->program_version eq "h3r") {
    print "ok 10\n";
}
else {
    print "not ok 10\n";
}

if ($analysis->gff_source eq "gff_test") {
    print "ok 11\n";
}
else {
    print "not ok 11\n";
}

if ($analysis->gff_feature eq "gff_feat") {
    print "ok 12\n";
}
else {
    print "not ok 12\n";
}

#Test the exist_analysis method
if ($feature_obj->exists_Analysis($analysis) == 1) {
    print "ok 13\n";
}
else {
    print "not ok 13";
}

#Set a feature object
my $out = Bio::EnsEMBL::FeatureFactory->new_feature_pair();
$out->set_all_fields(3,5,-1,18,'feat_temp','repeat',$contig->id,12,14,1,36,'feat_temp','repeat','h_temp');
$out->analysis($analysis);

#Call th  method to write this feature into the database
$feature_obj->write($contig,$out);

#Cal the method to get all similarities features in the contig
my @seq_feat = $contig->get_all_SimilarityFeatures;

if ((scalar(@seq_feat) != 2) && ($seq_feat[0]->id != 1)) {
print "not ok 14\n";
}
else {
    print "ok 14\n";
}

if ((scalar(@seq_feat) != 2) && ($seq_feat[1]->id != 2)) {
print "not ok 15\n";
}
else {
    print "ok 15\n";
}

@genome_hits = $feature_obj->find_GenomeHits("id_test");
if ((scalar(@genome_hits) != 1) && ($genome_hits[0]->start != 2)) {
    print "not ok 16\n";
}
else {
    print "ok 16\n";
}
#Delete all of the features from the contig 1
$feature_obj->delete(1);


$contig = $db->get_Contig("id_cont_test3");

#get again all of the similarityFeatures (there should not be there anymore...)
my @seq_feat1 = $contig->get_all_SimilarityFeatures;

if (scalar(@seq_feat1) != 0) {
print "not ok 17\n";
}
else {
    print "ok 17\n";
}


if ( ($seq_feat[0]->start == 3) && ($seq_feat[0]->end ==5) && ($seq_feat[0]->strand == -1)) {
    print "ok 18\n";
}
else {
    print "not ok 18\n";
    print STDERR $seq_feat[0]->start," ",$seq_feat[0]->end," ",$seq_feat[0]->strand,"\n";
}

my $anal2    = new Bio::EnsEMBL::Analysis(-db              => "feat2",
					  -db_version      => "3r",
					  -program         => "restgene",
					  -program_version => "5t",
					  -gff_source      => "gff_2",
					  -gff_feature     => "gff_3",
					  );
#write the analysis into the database.
$feature_obj->write_Analysis($anal2);

#Check if the the analysis has been written to the database
if ($feature_obj->exists_Analysis($anal2) == 2) {
    print "ok 19\n";
}
else {
    print "not ok 19\n";
}
exit;
