use lib 't';

BEGIN { $| = 1;
	use Test;
	plan tests => 2;
}


use MultiTestDB;
use TestUtils qw(test_getter_setter debug);

our $verbose = 0;

my $multi = MultiTestDB->new();
ok(1);

my $db = $multi->get_DBAdaptor('core');
my $pafa = $db->get_ProteinAlignFeatureAdaptor();

# list_dbIDs
debug("ProteinAlignFeatureAdaptor->list_dbIDs");
my $ids = $pafa->list_dbIDs();
ok (@{$ids});

