use strict;
use warnings;

use lib 't';

BEGIN { $| = 1;  
	use Test;
	plan tests => 7;
}


use TestUtils qw(debug test_getter_setter);
use Bio::EnsEMBL::KaryotypeBand;


#
#1 TEST - KaryotypeBand compiles
#
ok(1); 

#
# 2 test constructor
#
my $kb = Bio::EnsEMBL::KaryotypeBand->new;
ok($kb->isa('Bio::EnsEMBL::KaryotypeBand'));


#
# 3-7 test getter/setters
#
ok(test_getter_setter($kb, 'name', 'p.31'));
ok(test_getter_setter($kb, 'chr_name', 'X'));
ok(test_getter_setter($kb, 'start', 12_200_000));
ok(test_getter_setter($kb, 'end',   13_000_000));
ok(test_getter_setter($kb, 'stain', 'gpos50'));

