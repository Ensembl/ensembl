use strict;

use Bio::EnsEMBL::Test::TestUtils;

BEGIN { $| = 1;  
	use Test;
	plan tests => 10;
}

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::AssemblyExceptionFeature;

our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;

# get a core DBAdaptor
#
my $dba = $multi->get_DBAdaptor("core");
my $aefa = $dba->get_AssemblyExceptionFeatureAdaptor();
ok($aefa);

#
# 1 create a new AssemblyExceptionFeature
#
my $aef = new Bio::EnsEMBL::AssemblyExceptionFeature;
ok($aef);

#
# test the basic getter and setters
#

# start
ok(test_getter_setter($aef,'start',10));

# end
ok(test_getter_setter($aef,'end',14));

# type
ok(test_getter_setter($aef,'type', 'HAP'));

# check adaptor attaching
$aef->adaptor($aefa);
ok($aef->adaptor->isa('Bio::EnsEMBL::DBSQL::AssemblyExceptionFeatureAdaptor'));

# fetch all
my $chr_slice = $dba->get_SliceAdaptor->fetch_by_region('chromosome', 
                                                        '20_HAP1');
my @features = @{$aefa->fetch_all_by_Slice($chr_slice)};

ok(@features);
foreach my $f (@features) {
  debug( "Feature: " . $f->slice->seq_region_name . " " . 
         $f->start . " " . $f->end . " " . $f->type);
  my $as = $f->alternate_slice();
  debug(" Alternate slice: " . $as->seq_region_name . " " . 
        $as->start . " " . $as->end);
}

my ($f) = @features;
ok($f->display_id eq $f->alternate_slice->seq_region_name);


my $feat = $aefa->fetch_by_dbID(1);

ok($feat->dbID() == 1);

#check we can store assembly exception features
my $aef_store = new Bio::EnsEMBL::AssemblyExceptionFeature();
my $aef_store2 = new Bio::EnsEMBL::AssemblyExceptionFeature();
#get ref slice
my $ref_slice = $dba->get_SliceAdaptor->fetch_by_region('chromosome',20);
#prepare first object, the haplotype
$aef_store->start(1500);
$aef_store->end(35000);
$aef_store->type('HAP');
$aef_store->slice($chr_slice);
$aef_store->alternate_slice($ref_slice);
#prepare second object, the ref region to be substituted
$aef_store2->start(4500);
$aef_store2->end(38000);
$aef_store2->type('HAP REF');
$aef_store2->slice($ref_slice);
$aef_store2->alternate_slice($chr_slice);

my $asx_id = $aefa->store($aef_store,$aef_store2);

my $aef_new = $aefa->fetch_by_dbID($asx_id);

ok($aef_new->dbID == $asx_id);

$aefa->remove($aef_store);

