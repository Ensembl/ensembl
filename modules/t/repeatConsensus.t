use lib 't';
use TestUtils qw(test_getter_setter);

BEGIN { $| = 1;  
	use Test;
	plan tests => 13;
}

my $loaded = 0;
END {print "not ok 1\n" unless $loaded;}

my $verbose = 0;

use MultiTestDB;

my $multi = MultiTestDB->new();

$loaded = 1;

ok(1);

my $db = $multi->get_DBAdaptor( 'core' );

$repeat_c_ad = $db->get_RepeatConsensusAdaptor();

ok ($repeat_c_ad);

my $repeat_consensus = Bio::EnsEMBL::RepeatConsensus->new();

ok ($repeat_consensus);

ok (test_getter_setter($repeat_consensus,'length',10));
ok (test_getter_setter($repeat_consensus,'repeat_class','dummy'));
ok (test_getter_setter($repeat_consensus,'name','dummy'));
ok (test_getter_setter($repeat_consensus,'repeat_consensus','ATGCATGCAT'));
ok (test_getter_setter($repeat_consensus,'dbID',42));

ok ($repeat_consensus->desc eq 'class=dummy');
ok ($repeat_consensus->moltype eq 'dna');
ok ($repeat_consensus->alphabet eq 'dna');
ok ($repeat_consensus->seq eq 'ATGCATGCAT');

$repeat_c_ad->store($repeat_consensus);

ok(1);
