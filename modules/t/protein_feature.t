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
BEGIN { $| = 1; print "1..22\n"; 
	use vars qw($loaded); }
END {print "not ok 1\n" unless $loaded;}


use lib 't';
use EnsTestDB;
use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBSQL::Protein_Feature_Adaptor;

$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new();
    
# Load some data into the db
$ens_test->do_sql_file("t/protein_feature.dump");
    
# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;
print "ok 2\n";

my $protfeat = Bio::EnsEMBL::DBSQL::Protein_Feature_Adaptor->new($db);
my @transl;
my $id;
eval {
    @transl = $protfeat->Fetch_protfeature_by_translation("translation_id3");
};
if ($@) {
    print "not ok 3\n";
}
else {
    print "ok 3\n";
}

if ($transl[0]->feature1->start == 2) {
    print "ok 4\n";
}
else {
    print "not ok 4\n";
}

if ($transl[0]->feature1->end == 5) {
    print "ok 5\n";
}
else {
    print "not ok 5\n";
}

if ($transl[0]->feature1->score == 99.000) {
    print "ok 6\n";
}
else {
    print "not ok 6\n";
}

if ($transl[0]->feature1->analysis->id == 1) {
    print "ok 7\n";
}
else {
    print "not ok 7\n";
}

if ($transl[0]->feature1->seqname eq "translation_id3") {
    print "ok 8\n";
}
else {
    print "not ok 8\n";
}

if ($transl[0]->feature2->start == 4) {
    print "ok 9\n";
}
else {
    print "not ok 9\n";
}

if ($transl[0]->feature2->end == 7) {
    print "ok 10\n";
}
else {
    print "not ok 10\n";
}

if ($transl[0]->feature2->seqname eq "SP34") {
    print "ok 11\n";
}
else {
    print "not ok 11\n";
}

if (scalar(@transl) == 2) {
    print "ok 12\n";
}
else {
    print "not ok 12\n";
}

eval {
    $id = $protfeat->Fetch_protfeature_by_id(1);
};
if ($@) {
    print "not ok 13\n";
}
else {
    print "ok 13\n";
}

if ($id->feature1->start == $transl[0]->feature1->start) {
    print "ok 14\n";
}
else {
    print "not ok 14\n";
}

my $analysis = $protfeat->_feature_obj->get_Analysis(1);

my $feat1 = new Bio::EnsEMBL::SeqFeature ( -start => 1,                   
					   -end => 4,        
					   -score => 97,
					   -analysis => $analysis,
					   -seqname => "translation_id3");
   
my $feat2 = new Bio::EnsEMBL::SeqFeature (-start => 5,
					  -end => 6,
					  -analysis => $analysis,
					  -seqname => "SP408");
   

my $feature = new Bio::EnsEMBL::FeaturePair(-feature1 => $feat1,
					    -feature2 => $feat2);


eval {
$protfeat->write_Protein_feature($feature);
};

if ($@) {
    print "not ok 15\n";
}
else {
    print "ok 15\n";
}

@transl1 = $protfeat->Fetch_protfeature_by_translation("translation_id3");
if (scalar @transl1 == 3) {
    print "ok 16\n";
}
else {
    print "not ok 16\n";
}


eval {
$protfeat->delete("translation_id3");
};

if ($@) {
    print "not ok 17\n";
}
else {
    print "ok 17\n";
}




my $feat3 = new Bio::EnsEMBL::SeqFeature ( -start => 1,                   
					   -end => 4,        
					   -score => 97,
					   -analysis => $analysis,
					   -seqname => "translation_id3");
   
my $feat4 = new Bio::EnsEMBL::SeqFeature (-start => 5,
					  -end => 6,
					  -analysis => $analysis,
					  -seqname => "SP408");
   
my $feature1 = new Bio::EnsEMBL::FeaturePair(-feature1 => $feat3,
					     -feature2 => $feat4);


@transl1 = $protfeat->Fetch_protfeature_by_translation("translation_id3");

if (scalar @transl1 == 0) {
    print "ok 18\n";
}
else {
    print "not ok 18\n";
}

eval {
$protfeat->write_Protein_feature_by_translation("translation_id3",$feature,$feature1);
};

if ($@) {
    print "not ok 19\n";
}
else {
    print "ok 19\n";
}

my @transl2 = $protfeat->Fetch_protfeature_by_translation("translation_id3");

if (scalar @transl2 == 2) {
    print "ok 20\n";
}
else {
    print "not ok 20\n";
}

eval {
$protfeat->delete_by_id(1);
};

if ($@) {
    print "not ok 21\n";
}
else {
    print "ok 21\n";
}

eval {
    $id = $protfeat->Fetch_protfeature_by_id(1);
};
if ($@) {
    print "ok 22\n";
}
else {
    print "not ok 22\n";
}

