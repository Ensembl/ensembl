use lib 't';
use TestUtils qw(test_getter_setter);

BEGIN { $| = 1;  
	use Test;
	plan tests => 2;
}

use Bio::EnsEMBL::MarkerFeature;

#
# 1 create a new Simplefeature
#
$mf = new Bio::EnsEMBL::MarkerFeature;
ok($mf);


#
# 2 test the only method getter and setters
#

# 2 start
my $test = 0;
$test = test_getter_setter($mf,'marker_name','dummy_marker');
ok($test);
