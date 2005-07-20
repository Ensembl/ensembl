use strict;

use Bio::EnsEMBL::Test::TestUtils;

BEGIN { $| = 1;
	use Test;
	plan tests => 9;
}

use Bio::EnsEMBL::Test::MultiTestDB;
use Bio::EnsEMBL::RegulatoryFeature;
use Bio::EnsEMBL::RegulatoryFactor;


our $verbose = 0;

my $multi = Bio::EnsEMBL::Test::MultiTestDB->new;

# get a core DBAdaptor
#
my $dba = $multi->get_DBAdaptor("core");
my $rfa = $dba->get_RegulatoryFeatureAdaptor;

#
# 1 create a new RegulatoryFeature
#
my $rf = new Bio::EnsEMBL::RegulatoryFeature;
ok($rf);


#
# 2-7 test the basic getter and setters
#

# 2 start
ok(test_getter_setter($rf,'start',10));

# 3 end
ok(test_getter_setter($rf,'end',14));

# 4 strand
ok(test_getter_setter($rf,'strand',1));

# 5 factor
#my $rm = Bio::EnsEMBL::RegulatoryFactor->new(-NAME => 'Joe',
#					    -TYPE => 'promoter');
#ok(test_getter_setter($rf,'factor',$rm));

# 6 influence
#ok(test_getter_setter($rf,'influence','positive'));

# 7 dbID
ok(test_getter_setter($rf,'dbID',20));

#
# 8 check adaptor attaching
#
$rf->adaptor($rfa);
ok($rf->adaptor->isa('Bio::EnsEMBL::DBSQL::RegulatoryFeatureAdaptor'));

# 9 check retrieving regulated transcripts
my $rf1 = $rfa->fetch_by_dbID(1);
ok(@{$rf1->regulated_transcripts()} == 1);

# 10 check get_transcripts_regulated_by_same_factor
ok(@{$rf1->transcripts_regulated_by_same_factor()} == 1);



