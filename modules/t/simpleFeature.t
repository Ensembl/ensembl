use strict;

use lib 't';
use TestUtils qw(test_getter_setter debug);

BEGIN { $| = 1;  
	use Test;
	plan tests => 10;
}

use MultiTestDB;
use Bio::EnsEMBL::SimpleFeature;


our $verbose = 1;

my $multi = MultiTestDB->new;
 
# get a core DBAdaptor
#
my $dba = $multi->get_DBAdaptor("core");
my $sfa = $dba->get_SimpleFeatureAdaptor;

#
# 1 create a new Simplefeature
#
my $sf = new Bio::EnsEMBL::SimpleFeature;
ok($sf);


#
# 2-7 test the basic getter and setters
#

# 2 start
ok(test_getter_setter($sf,'start',10));

# 3 end
ok(test_getter_setter($sf,'end',14));

# 4 strand
ok(test_getter_setter($sf,'strand',1));

# 5 score
ok(test_getter_setter($sf,'score',42));

# 6 display_label
ok(test_getter_setter($sf,'display_label','dummy_label'));

# 7 dbID
ok(test_getter_setter($sf,'dbID',42));



#
# 9 check adaptor attaching
#
$sf->adaptor($sfa);
ok($sf->adaptor->isa('Bio::EnsEMBL::DBSQL::SimpleFeatureAdaptor'));


my $slice = $dba->get_SliceAdaptor->fetch_by_region('chromosome', '20');

my $features = $sfa->fetch_all_by_Slice($slice);

foreach my $feature (@$features) {
  debug("\n\nfeature start = " . $feature->start());
  debug("feature end = " . $feature->end());
  debug("feature strand = " . $feature->strand());
  debug("feature display_label = " . $feature->display_label());
  debug("feature score = " . $feature->score());
}

# List_dbidx
my $ids = $sfa->list_dbIDs();
ok (@{$ids});
