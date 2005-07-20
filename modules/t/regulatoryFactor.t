use strict;
use Bio::EnsEMBL::Test::TestUtils;

use Bio::EnsEMBL::RegulatoryMotif;

BEGIN { $| = 1;
	use Test;
	plan tests => 7;
}

my $verbose = 0;

#
# Test constructor
#
my $rm = Bio::EnsEMBL::RegulatoryMotif->new(-NAME => 'Joe',
					    -TYPE => 'promoter');

ok($rm);
ok(ref($rm));
ok($rm->isa('Bio::EnsEMBL::RegulatoryMotif'));

ok($rm->name eq 'Joe');
ok($rm->type eq 'promoter');

ok(test_getter_setter($rm,'name','Fred'));
ok(test_getter_setter($rm,'type','miRNA_target'));
