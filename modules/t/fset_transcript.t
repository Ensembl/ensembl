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


## We start with some black magic to print on failure.
BEGIN { $| = 1; print "1..9\n"; 
	use vars qw($loaded); }

END {print "not ok 1\n" unless $loaded;}

use Bio::EnsEMBL::DBSQL::Obj;
use Bio::EnsEMBL::DBLoader;
use lib 't';
use EnsTestDB;
$loaded = 1;
print "ok 1\n";    # 1st test passes.
    
my $ens_test = EnsTestDB->new();
    
# Load some data into the db
$ens_test->do_sql_file("t/donor.dump");


# Get an EnsEMBL db object for the test db
my $db = $ens_test->get_DBSQL_Obj;
print "ok 2\n";

#Get a new feature_obj object
my $feature_obj=Bio::EnsEMBL::DBSQL::Feature_Obj->new($db);
print "ok 3\n";



if ($feature_obj->get_PredictionFeature_by_id("41")->isa ("Bio::EnsEMBL::SeqFeatureI"))
{print "ok 4\n";}
else {
    print "not ok 4\n";
}

eval {
$feature_obj->get_PredictionFeature_by_id("wrong_id");
};

if ($@){print "ok 5\n";}
else { print "not ok 5\n";}

if ($feature_obj->get_PredictionFeature_as_Transcript(40)->isa ("Bio::EnsEMBL::Transcript"))
{print  "ok 6\n";}
else { print "not ok 6\n";}


eval {
$feature_obj->get_PredictionFeature_as_Transcript("wrong_id");
};

if ($@){print "ok 7\n";}
else { print "not ok 7\n";}



if ($db->get_PredictionFeature_as_Transcript(40)->isa ("Bio::EnsEMBL::Transcript"))
{print  "ok 8\n";}
else { print "not ok 8\n";}


eval {
$db->get_PredictionFeature_as_Transcript("wrong_id");
};

if ($@){print "ok 9\n";}
else { print "not ok 9\n";}














