use strict;

use lib 't';
use TestUtils qw(test_getter_setter debug);

BEGIN { $| = 1;  
	use Test;
	plan tests => 8;
}

use MultiTestDB;
use Bio::EnsEMBL::AssemblyExceptionFeature;


our $verbose = 0;

my $multi = MultiTestDB->new;

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

my $f = (@features);
ok($f->display_id eq $f->alternate_slice->seq_region_name);
